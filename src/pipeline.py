from pathlib import Path
from typing import Dict

from src.data.gdc_api import fetch_cases
from src.cohort.build_survival import build_survival_table
from src.cohort.load_curated_survival import load_tcga_cdr_survival

def run_pipeline(cfg: Dict, run_dir: Path) -> None:
    logs_dir = run_dir / "logs"
    logs_dir.mkdir(exist_ok=True)

    project = cfg["data"]["cancer"]
    processed_dir = Path(cfg["data"]["processed_dir"])
    processed_dir.mkdir(parents=True, exist_ok=True)

    fields = [
        "submitter_id",
        "demographic.vital_status",
        "demographic.days_to_death",
        "demographic.days_to_last_follow_up",
    ]

    (logs_dir / "pipeline.log").write_text("Fetching clinical cases from GDC...\n")

    cases = fetch_cases(project=project, fields=fields)
    df_surv, df_flags = build_survival_table(cases)

    surv_path = processed_dir / "survival_table.csv"
    flags_path = processed_dir / "survival_flags.csv"
    df_surv.to_csv(surv_path, index=False)
    df_flags.to_csv(flags_path, index=False)

    (logs_dir / "pipeline.log").write_text(
        (logs_dir / "pipeline.log").read_text()
        + f"Wrote {len(df_surv)} usable patients to {surv_path}\n"
        + f"Wrote flags to {flags_path}\n"
    )

        # --- Curated survival (TCGA-CDR) ---
    curated_path = Path("data/curated/tcga_cdr_survival.xlsx")
    endpoint = cfg["data"].get("endpoint", "OS")
    cancer_type = cfg["data"].get("cancer_short", "LUAD")

    df_curated = load_tcga_cdr_survival(
        xlsx_path=curated_path,
        cancer_type=cancer_type,
        endpoint=endpoint,
    )

    processed_dir = Path(cfg["data"]["processed_dir"])
    processed_dir.mkdir(parents=True, exist_ok=True)

    curated_out = processed_dir / "survival_table_curated.csv"
    df_curated.to_csv(curated_out, index=False)

    (logs_dir / "pipeline.log").write_text(
        (logs_dir / "pipeline.log").read_text()
        + f"Curated survival: wrote {len(df_curated)} rows to {curated_out}\n"
        + f"Curated event counts:\n{df_curated['event'].value_counts().to_string()}\n"
    )