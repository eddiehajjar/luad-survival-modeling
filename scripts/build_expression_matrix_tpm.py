from __future__ import annotations

from pathlib import Path
import pandas as pd


RAW_DIR = Path("data/raw/expression/HTSeq-FPKM")
OUT_PATH = Path("data/processed/expression_tpm_matrix.csv")


def read_one(path: Path) -> pd.Series:
    # Skip comment lines starting with '#'
    df = pd.read_csv(
        path,
        sep=r"\s+",
        comment="#",
        engine="python",
        dtype={"gene_id": str},
    )

    # Keep only real genes
    df = df[df["gene_id"].str.startswith("ENSG")].copy()

    # Use TPM
    s = df.set_index("gene_id")["tpm_unstranded"]
    s.name = path.stem  # file_id
    return s


def main():
    files = sorted(RAW_DIR.glob("*.tsv"))
    if not files:
        raise SystemExit(f"No .tsv files found in {RAW_DIR}")

    print("Found expression files:", len(files))

    series_list = []
    for p in files:
        try:
            series_list.append(read_one(p))
        except Exception as e:
            print("Failed:", p.name, "error:", e)

    mat = pd.concat(series_list, axis=1).T
    mat.index.name = "file_id"

    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    mat.to_csv(OUT_PATH)

    print("Wrote TPM matrix to:", OUT_PATH)
    print("Matrix shape:", mat.shape)


if __name__ == "__main__":
    main()