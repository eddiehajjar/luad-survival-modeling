from __future__ import annotations

from pathlib import Path
import pandas as pd

from src.data.gdc_expression import list_luad_expression_files


EXPR_PATH = Path("data/processed/expression_tpm_matrix.csv")
SURV_PATH = Path("data/processed/survival_table_curated.csv")
OUT_PATH = Path("data/processed/luad_survival_expression.csv")


def main():
    # Load expression matrix
    expr = pd.read_csv(EXPR_PATH)
    print("Expression matrix:", expr.shape)

    # Build file_id -> patient_id manifest
    files = list_luad_expression_files(workflow_type=None, sample_type="Primary Tumor", max_files=5000)
    manifest = pd.DataFrame(
        [{"file_id": f.file_id, "patient_id": f.patient_id} for f in files]
    )

    print("Manifest rows:", len(manifest))

    # Merge expression with manifest
    expr = expr.merge(manifest, on="file_id", how="inner")
    print("After adding patient_id:", expr.shape)

    # --- Deduplicate patients: keep ONE tumor profile per patient ---
    # Rule: keep the row with the largest total TPM across all genes
    gene_cols = [c for c in expr.columns if c.startswith("ENSG")]
    expr["tpm_sum"] = expr[gene_cols].sum(axis=1)

    expr = (
        expr.sort_values(["patient_id", "tpm_sum"], ascending=[True, False])
            .drop_duplicates(subset=["patient_id"], keep="first")
            .drop(columns=["tpm_sum"])
    )

    print("After deduping patients:", expr.shape)

    # Load survival table
    surv = pd.read_csv(SURV_PATH)
    print("Survival table:", surv.shape)

    # Merge expression with survival
    final = expr.merge(surv, on="patient_id", how="inner")

    # Drop file_id (not needed downstream)
    final = final.drop(columns=["file_id"])

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    final.to_csv(OUT_PATH, index=False)

    print("Final aligned table written to:", OUT_PATH)
    print("Final shape:", final.shape)


if __name__ == "__main__":
    main()