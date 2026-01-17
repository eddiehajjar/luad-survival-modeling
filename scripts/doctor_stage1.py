from __future__ import annotations

from pathlib import Path
import sys
import pandas as pd


def ok(msg: str) -> None:
    print(f"[OK]   {msg}")


def warn(msg: str) -> None:
    print(f"[MISS] {msg}")


def main() -> int:
    root = Path(".")
    problems = 0

    # --- Survival labels (curated) ---
    surv_path = root / "data/processed/survival_table_curated.csv"
    if surv_path.exists():
        df = pd.read_csv(surv_path)
        ok(f"Found curated survival table: {surv_path} (rows={len(df)})")
        for col in ["patient_id", "time_days", "event"]:
            if col not in df.columns:
                warn(f"Curated survival missing column: {col}")
                problems += 1
        if "event" in df.columns:
            vc = df["event"].value_counts(dropna=False).to_dict()
            ok(f"Event counts: {vc}")
    else:
        warn(f"Missing curated survival table: {surv_path}")
        problems += 1

    # --- Raw expression files inventory ---
    expr_dir_a = root / "data/raw/expression/HTSeq-FPKM"
    expr_dir_b = root / "data/raw/expression/STAR-Counts"

    expr_dir = None
    if expr_dir_a.exists():
        expr_dir = expr_dir_a
    elif expr_dir_b.exists():
        expr_dir = expr_dir_b

    if expr_dir is None:
        warn("Missing expression directory: data/raw/expression/HTSeq-FPKM or data/raw/expression/STAR-Counts")
        problems += 1
        expr_files = []
    else:
        expr_files = sorted(expr_dir.glob("*.tsv"))
        ok(f"Expression directory: {expr_dir} (tsv_files={len(expr_files)})")
        if len(expr_files) == 0:
            warn("No expression .tsv files found (need to download LUAD expression files)")
            problems += 1

    # --- Expression matrix (optional, next milestone) ---
    mat_path = root / "data/processed/expression_tpm_matrix.csv"
    if mat_path.exists():
        # don’t read full file if huge; just read header and a few rows
        dfm = pd.read_csv(mat_path, nrows=5)
        ok(f"Found TPM matrix: {mat_path} (preview rows=5, cols={len(dfm.columns)})")
    else:
        warn(f"Missing TPM matrix: {mat_path} (will be created after full download)")
        # not counting as a hard fail yet; depends on where you are

    # --- Final aligned modeling table (later milestone) ---
    aligned_path = root / "data/processed/luad_survival_expression.csv"
    if aligned_path.exists():
        dfa = pd.read_csv(aligned_path, nrows=5)
        ok(f"Found aligned modeling table: {aligned_path} (preview cols={len(dfa.columns)})")
    else:
        warn(f"Missing aligned modeling table: {aligned_path} (created after mapping file_id→patient_id and joining)")
        # not counting as hard fail yet

    print("\n--- Summary ---")
    if problems == 0:
        ok("All required Stage 1 prerequisites found (at your current step).")
        return 0
    else:
        warn(f"{problems} required items missing or incomplete.")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())