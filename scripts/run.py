from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from src.pipeline import run_pipeline
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

import yaml


def load_config(path: str) -> dict:
    with open(path, "r") as f:
        return yaml.safe_load(f)


def make_run_dir(runs_dir: str, run_name: str) -> Path:
    ts = datetime.now().strftime("%Y-%m-%d_%H%M%S")
    run_dir = Path(runs_dir) / f"{ts}_{run_name}"
    run_dir.mkdir(parents=True, exist_ok=False)
    return run_dir


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    args = parser.parse_args()

    cfg = load_config(args.config)
    run_name = cfg["project"]["run_name"]
    runs_dir = cfg["outputs"]["runs_dir"]

    run_dir = make_run_dir(runs_dir, run_name)

    # Save an exact copy of the config used for this run
    (run_dir / "config_used.yaml").write_text(yaml.safe_dump(cfg, sort_keys=False))

    # Save run metadata
    meta = {
        "created_at": datetime.now().isoformat(timespec="seconds"),
        "config_path": args.config,
        "run_dir": str(run_dir),
    }
    (run_dir / "run_meta.json").write_text(json.dumps(meta, indent=2))
    run_pipeline(cfg, run_dir)

    print(f"Run directory created: {run_dir}")
    print(f"Run completed. See logs at: {run_dir / 'logs' / 'pipeline.log'}")

if __name__ == "__main__":
    main()