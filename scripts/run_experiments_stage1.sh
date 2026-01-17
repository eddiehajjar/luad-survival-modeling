#!/usr/bin/env bash
set -euo pipefail

# Always run from repo root
cd "$(dirname "$0")/.."

TS="$(date +%Y-%m-%d_%H%M%S)"
OUTROOT="artifacts/experiments/$TS"
LOGDIR="$OUTROOT/logs"
mkdir -p "$LOGDIR"

echo "[INFO] Writing results under: $OUTROOT"
echo "[INFO] Logs under: $LOGDIR"

# Helper: run one experiment and tee logs
run_one () {
  local NAME="$1"
  shift
  local LOGFILE="$LOGDIR/${NAME}.log"

  echo "==================================================" | tee -a "$LOGFILE"
  echo "[RUN] $NAME" | tee -a "$LOGFILE"
  echo "[CMD] PYTHONPATH=. python scripts/train_cox_cv.py $*" | tee -a "$LOGFILE"
  echo "[TIME] $(date)" | tee -a "$LOGFILE"
  echo "==================================================" | tee -a "$LOGFILE"

  PYTHONPATH=. python scripts/train_cox_cv.py \
    --outdir "$OUTROOT/$NAME" \
    "$@" 2>&1 | tee -a "$LOGFILE"

  echo "[DONE] $NAME at $(date)" | tee -a "$LOGFILE"
}

# --- Experiments ---
# 1) Your current run setup (for comparison)
run_one "cox_5000_pen0p1_f5" --top_genes 5000 --folds 5 --penalizer 0.1 --l1_ratio 0.0

# 2) Stronger regularization (often improves generalization)
run_one "cox_5000_pen1_f5"   --top_genes 5000 --folds 5 --penalizer 1.0 --l1_ratio 0.0
run_one "cox_5000_pen5_f5"   --top_genes 5000 --folds 5 --penalizer 5.0 --l1_ratio 0.0

# 3) Fewer genes (often improves generalization at this sample size)
run_one "cox_3000_pen1_f5"   --top_genes 3000 --folds 5 --penalizer 1.0 --l1_ratio 0.0
run_one "cox_2000_pen1_f5"   --top_genes 2000 --folds 5 --penalizer 1.0 --l1_ratio 0.0

echo "[ALL DONE] Completed at $(date)"
echo "[RESULTS] $OUTROOT"