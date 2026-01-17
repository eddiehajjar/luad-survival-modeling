from __future__ import annotations
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pandas as pd

def _get(d: Dict[str, Any], path: str):
    cur = d
    for part in path.split("."):
        if cur is None:
            return None
        if isinstance(cur, dict):
            cur = cur.get(part)
        else:
            return None
    return cur

def build_survival_table(cases: List[Dict[str, Any]]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    flags = []

    for c in cases:
        pid = c.get("submitter_id") or c.get("case_id")

        vital = _get(c, "demographic.vital_status")
        dtd = _get(c, "demographic.days_to_death")
        dtlfu = _get(c, "demographic.days_to_last_follow_up")

        # Normalize numeric fields
        def to_float(x):
            try:
                return float(x) if x is not None else None
            except Exception:
                return None

        dtd = to_float(dtd)
        dtlfu = to_float(dtlfu)

        # Core rule: prefer observed death time if present.
        # event=1 if days_to_death exists and is positive
        event = 1 if (dtd is not None and dtd > 0) else 0
        time_days = dtd if event == 1 else dtlfu

        # Flag weird/inconsistent cases for later inspection (do NOT hide them)
        f = {"patient_id": pid, "vital_status": vital, "days_to_death": dtd, "days_to_last_follow_up": dtlfu}
        if time_days is None:
            f["flag"] = "missing_time"
        elif time_days <= 0:
            f["flag"] = "nonpositive_time"
        elif vital is not None and str(vital).lower() == "dead" and event == 0:
            f["flag"] = "dead_but_no_days_to_death"
        elif vital is not None and str(vital).lower() == "alive" and event == 1:
            f["flag"] = "alive_but_has_days_to_death"
        else:
            f["flag"] = ""
        flags.append(f)

        rows.append({"patient_id": pid, "time_days": time_days, "event": event})

    df = pd.DataFrame(rows)
    flag_df = pd.DataFrame(flags)

    # Exclude unusable rows for modeling, but keep a record of what was excluded
    usable = df["time_days"].notna() & (df["time_days"] > 0)
    df_clean = df.loc[usable].copy()

    return df_clean, flag_df