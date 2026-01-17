from __future__ import annotations

from pathlib import Path
import pandas as pd


def load_tcga_cdr_survival(
    xlsx_path: Path,
    cancer_type: str = "LUAD",
    endpoint: str = "OS",
) -> pd.DataFrame:
    """
    Load curated survival labels from the TCGA Clinical Data Resource (TCGA-CDR).

    Returns a DataFrame with columns:
      patient_id, time_days, event
    """
    df = pd.read_excel(xlsx_path, sheet_name="TCGA-CDR")

    time_col = f"{endpoint}.time"
    event_col = endpoint

    # Filter to cancer type (e.g., LUAD)
    df = df[df["type"] == cancer_type].copy()

    # Standardize schema
    out = df.rename(
        columns={
            "bcr_patient_barcode": "patient_id",
            time_col: "time_days",
            event_col: "event",
        }
    )[["patient_id", "time_days", "event"]]

    # Clean
    out = out.dropna(subset=["time_days", "event"])
    out = out[out["time_days"] > 0].copy()
    out["event"] = out["event"].astype(int)

    return out