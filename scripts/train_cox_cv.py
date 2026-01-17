from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.statistics import logrank_test
from sklearn.model_selection import StratifiedKFold

import matplotlib.pyplot as plt


def log(msg: str) -> None:
    print(msg, flush=True)


def load_data(path: str) -> pd.DataFrame:
    # Parquet is preferred; CSV fallback if you want
    if path.endswith(".parquet"):
        return pd.read_parquet(path)
    return pd.read_csv(path)


def variance_filter_train_only(X_train: np.ndarray, gene_cols: list[str], top_k: int) -> list[str]:
    # variance computed ONLY on training fold -> no leakage
    variances = X_train.var(axis=0)
    top_idx = np.argsort(variances)[::-1][:top_k]
    return [gene_cols[i] for i in top_idx]


def fit_cox(train_df: pd.DataFrame, penalizer: float, l1_ratio: float) -> CoxPHFitter:
    cph = CoxPHFitter(penalizer=penalizer, l1_ratio=l1_ratio)
    cph.fit(train_df, duration_col="time_days", event_col="event")
    return cph


def c_index_on(cph: CoxPHFitter, df: pd.DataFrame) -> float:
    # lifelines score gives concordance_index for CoxPHFitter
    return float(cph.score(df, scoring_method="concordance_index"))


def km_plot(df: pd.DataFrame, risk: pd.Series, out_path: Path) -> dict:
    # High vs Low by median risk score
    median = float(risk.median())
    high = df[risk >= median]
    low = df[risk < median]

    kmf = KaplanMeierFitter()

    plt.figure()
    kmf.fit(low["time_days"], event_observed=low["event"], label="Low risk")
    ax = kmf.plot()

    kmf.fit(high["time_days"], event_observed=high["event"], label="High risk")
    kmf.plot(ax=ax)

    plt.xlabel("Days")
    plt.ylabel("Survival probability")
    plt.title("Kaplan–Meier: High vs Low risk (median split)")
    plt.tight_layout()
    plt.savefig(out_path)
    plt.close()

    # Log-rank test
    lr = logrank_test(
        low["time_days"], high["time_days"],
        event_observed_A=low["event"], event_observed_B=high["event"]
    )
    return {
        "median_risk_threshold": median,
        "n_low": int(len(low)),
        "n_high": int(len(high)),
        "logrank_p_value": float(lr.p_value),
        "test_statistic": float(lr.test_statistic),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data", default="data/processed/luad_survival_expression.parquet")
    parser.add_argument("--outdir", default="artifacts/model_runs/cox_cv")
    parser.add_argument("--top_genes", type=int, default=5000)
    parser.add_argument("--penalizer", type=float, default=0.1)
    parser.add_argument("--l1_ratio", type=float, default=0.0)  # 0=ridge, 1=lasso
    parser.add_argument("--folds", type=int, default=5)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    t0 = time.time()
    outdir = Path(args.outdir) / time.strftime("%Y-%m-%d_%H%M%S")
    outdir.mkdir(parents=True, exist_ok=False)
    (outdir / "plots").mkdir(parents=True, exist_ok=True)
    (outdir / "tables").mkdir(parents=True, exist_ok=True)

    log(f"[INFO] Output dir: {outdir}")

    df = load_data(args.data)
    log(f"[INFO] Loaded data: {df.shape}")

    # Basic sanity checks
    required = {"patient_id", "time_days", "event"}
    missing = required - set(df.columns)
    if missing:
        raise SystemExit(f"Missing required columns: {missing}")

    if df["patient_id"].duplicated().any():
        raise SystemExit("Duplicate patient_id detected. Deduplicate before modeling.")

    # Identify gene columns
    gene_cols = [c for c in df.columns if c.startswith("ENSG")]
    if len(gene_cols) < args.top_genes:
        raise SystemExit(f"Requested top_genes={args.top_genes} but only {len(gene_cols)} gene columns found.")

    # y
    y = df[["time_days", "event"]].copy()

    # X: log1p(TPM) as float32 for speed/memory
    log("[INFO] Building log1p expression matrix (float32)...")
    X = np.log1p(df[gene_cols].to_numpy(dtype=np.float32))
    log(f"[INFO] X shape: {X.shape}")

    # CV splits stratified by event indicator (keeps event ratio similar per fold)
    skf = StratifiedKFold(n_splits=args.folds, shuffle=True, random_state=args.seed)

    fold_rows = []
    coef_paths = []
    gene_list_paths = []

    log(f"[INFO] Starting {args.folds}-fold CV (top_genes selected inside each fold)...")

    for fold, (tr_idx, va_idx) in enumerate(skf.split(X, y["event"]), start=1):
        fold_t0 = time.time()

        X_tr, X_va = X[tr_idx], X[va_idx]
        y_tr, y_va = y.iloc[tr_idx].reset_index(drop=True), y.iloc[va_idx].reset_index(drop=True)

        # Feature selection on training only (no leakage)
        top_genes = variance_filter_train_only(X_tr, gene_cols, args.top_genes)

        # Build train/val dataframes for lifelines
        tr_df = pd.concat(
            [y_tr, pd.DataFrame(X_tr[:, [gene_cols.index(g) for g in top_genes]], columns=top_genes)],
            axis=1
        )
        va_df = pd.concat(
            [y_va, pd.DataFrame(X_va[:, [gene_cols.index(g) for g in top_genes]], columns=top_genes)],
            axis=1
        )

        # Fit
        cph = fit_cox(tr_df, penalizer=args.penalizer, l1_ratio=args.l1_ratio)

        # Scores
        train_c = float(cph.concordance_index_)  # training concordance
        val_c = c_index_on(cph, va_df)

        # Save artifacts per fold
        coef_path = outdir / "tables" / f"fold_{fold:02d}_coefs.csv"
        gene_path = outdir / "tables" / f"fold_{fold:02d}_genes.txt"

        cph.params_.to_csv(coef_path, header=["coef"])
        Path(gene_path).write_text("\n".join(top_genes))

        coef_paths.append(str(coef_path))
        gene_list_paths.append(str(gene_path))

        fold_rows.append({
            "fold": fold,
            "train_n": int(len(tr_df)),
            "val_n": int(len(va_df)),
            "train_c_index": train_c,
            "val_c_index": val_c,
            "runtime_s": float(time.time() - fold_t0),
        })

        log(f"[FOLD {fold}] train_c={train_c:.4f} val_c={val_c:.4f}  (t={time.time()-fold_t0:.1f}s)")

    metrics = pd.DataFrame(fold_rows)
    metrics_path = outdir / "tables" / "cv_metrics.csv"
    metrics.to_csv(metrics_path, index=False)

    summary = {
        "data": args.data,
        "top_genes": args.top_genes,
        "penalizer": args.penalizer,
        "l1_ratio": args.l1_ratio,
        "folds": args.folds,
        "seed": args.seed,
        "val_c_index_mean": float(metrics["val_c_index"].mean()),
        "val_c_index_std": float(metrics["val_c_index"].std(ddof=1)),
        "train_c_index_mean": float(metrics["train_c_index"].mean()),
        "train_c_index_std": float(metrics["train_c_index"].std(ddof=1)),
        "coef_paths": coef_paths,
        "gene_list_paths": gene_list_paths,
        "total_runtime_s": float(time.time() - t0),
    }

    (outdir / "run_summary.json").write_text(json.dumps(summary, indent=2))

    log("\n--- CV Summary ---")
    log(f"Val C-index mean: {summary['val_c_index_mean']:.4f}")
    log(f"Val C-index std : {summary['val_c_index_std']:.4f}")
    log(f"Train C-index mean: {summary['train_c_index_mean']:.4f}")
    log(f"Train C-index std : {summary['train_c_index_std']:.4f}")
    log(f"[INFO] Saved metrics: {metrics_path}")
    log(f"[INFO] Saved summary: {outdir / 'run_summary.json'}")

    # Fit one final model on ALL data (for interpretation + KM plot)
    # Feature selection on full data is fine here because this model is not being “evaluated”; it’s for interpretation.
    log("\n[INFO] Fitting final model on full dataset for interpretation...")
    top_genes_full = variance_filter_train_only(X, gene_cols, args.top_genes)
    full_df = pd.concat(
        [y, pd.DataFrame(X[:, [gene_cols.index(g) for g in top_genes_full]], columns=top_genes_full)],
        axis=1
    )
    cph_final = fit_cox(full_df, penalizer=args.penalizer, l1_ratio=args.l1_ratio)

    final_coef_path = outdir / "tables" / "final_model_coefs.csv"
    cph_final.params_.to_csv(final_coef_path, header=["coef"])
    (outdir / "tables" / "final_model_genes.txt").write_text("\n".join(top_genes_full))

    # Risk scores for KM plot
    risk_scores = cph_final.predict_partial_hazard(full_df).rename("risk_score")

    km_info = km_plot(
        df=y.join(risk_scores),
        risk=risk_scores,
        out_path=outdir / "plots" / "km_high_vs_low.png",
    )

    (outdir / "km_info.json").write_text(json.dumps(km_info, indent=2))

    log(f"[INFO] Saved final coefs: {final_coef_path}")
    log(f"[INFO] Saved KM plot  : {outdir / 'plots' / 'km_high_vs_low.png'}")
    log(f"[INFO] Log-rank p-value: {km_info['logrank_p_value']:.3e}")
    log(f"[DONE] Total runtime: {summary['total_runtime_s']:.1f}s")


if __name__ == "__main__":
    main()