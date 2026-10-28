#!/bin/bash
# Unified MVMR runner for primary / gscan / extended phases.
#
# Usage:
#   cd /home/hello/influenza_mental_genomics
#   bash scripts/MVMR.run.sh <phase> [--dry-run]
#
# Phases:
#   primary   : clinical-diagnosis SUD MVMR (Tobacco/Alcohol/Opiate on influenza)
#   gscan     : GSCAN smoking-behaviour MVMR
#   extended  : extended infection-outcome MVMR
#
# Optional dry-run (print commands without executing):
#   bash scripts/MVMR.run.sh primary --dry-run

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"

PHASE="${1:-}"
if [[ -z "$PHASE" || "$PHASE" == "--dry-run" ]]; then
  echo "Usage: $0 <phase> [--dry-run]"
  echo ""
  echo "Available phases:"
  echo "  primary   : clinical-diagnosis SUD MVMR"
  echo "  gscan     : GSCAN smoking-behaviour MVMR"
  echo "  extended  : extended infection-outcome MVMR"
  exit 1
fi

DRY_RUN=0
if [[ "${2:-}" == "--dry-run" ]]; then
  DRY_RUN=1
fi

CONFIG="$SCRIPT_DIR/MVMR.${PHASE}.config.csv"
LOG_DIR="${MVMR_LOG_DIR:-$SCRIPT_DIR/logs}"
RSCRIPT="/tools/R-4.6.1/bin/Rscript"
SCRIPT="$SCRIPT_DIR/MVMR.single.R"

if [[ ! -f "$CONFIG" ]]; then
  echo "ERROR: config file not found: $CONFIG"
  exit 1
fi

mkdir -p "$LOG_DIR"

# Common environment (allow overrides via exported env vars)
export MVMR_PHASE="$PHASE"
export MVMR_INFO_FILE=${MVMR_INFO_FILE:-/home/hello/influenza_mental_genomics/GWAS.infos.xlsx}
export MVMR_DATA_DIR=${MVMR_DATA_DIR:-/home/hello/influenza_mental_genomics/data}
export MVMR_OUT_DIR=${MVMR_OUT_DIR:-/home/hello/influenza_mental_genomics/results/MVMR/${PHASE}}
export MVMR_PLINK_BIN=${MVMR_PLINK_BIN:-/tools/plink-1.90/plink}
export MVMR_BFILE_EUR=${MVMR_BFILE_EUR:-/tools/plink-1.90/1kg.v3/EUR}
export MVMR_HAPMAP3_SNPS=${MVMR_HAPMAP3_SNPS:-}

mkdir -p "$MVMR_OUT_DIR"

echo "============================================================"
echo "MVMR unified runner"
echo "============================================================"
echo "Phase:        $PHASE"
echo "Config:       $CONFIG"
echo "Data dir:     $MVMR_DATA_DIR"
echo "Output dir:   $MVMR_OUT_DIR"
echo "LD reference: $MVMR_BFILE_EUR"
echo "Dry-run:      $DRY_RUN"
echo "============================================================"

# Count models
n_models=$(tail -n +2 "$CONFIG" | awk -F',' 'NF>=2 && $1!="" {print}' | wc -l)
echo "Models to process: $n_models"
echo ""

# Concurrency limit (default 6; override with MVMR_MAX_JOBS)
MAX_JOBS=${MVMR_MAX_JOBS:-6}
running_jobs=0

# Skip header and process each row
# Use process substitution so variables are updated in the current shell
while IFS=, read -r model_name panel outcome_id exposure_ids covariate_ids p_thresh clump_kb clump_r2 note; do
  # Skip empty lines
  [[ -z "$model_name" ]] && continue

  # Config uses ';' inside covariate_ids, but MVMR.single.R expects comma-separated env var
  covariate_ids_comma="${covariate_ids//;/,}"

  export MVMR_MODEL_NAME="$model_name"
  export MVMR_OUTCOME_ID="$outcome_id"
  export MVMR_EXPOSURE_IDS="$exposure_ids"
  export MVMR_COVARIATE_IDS="$covariate_ids_comma"
  export MVMR_PANEL="$panel"
  export MVMR_P_THRESH="$p_thresh"
  export MVMR_CLUMP_KB="$clump_kb"
  export MVMR_CLUMP_R2="$clump_r2"

  LOG_FILE="$LOG_DIR/MVMR.${PHASE}.${model_name}.log"
  OUT_DIR="$MVMR_OUT_DIR/$model_name"

  # ===== Skip logic: check if results already exist =====
  if ls "$OUT_DIR"/*MVMR_IVW_results.csv 1> /dev/null 2>&1; then
    echo "[skip] $model_name (results already exist in $OUT_DIR)"
    continue
  fi
  # ======================================================

  if [[ $DRY_RUN -eq 1 ]]; then
    echo "[dry-run] Would submit: $model_name (panel=$panel, p=$p_thresh, kb=$clump_kb, r2=$clump_r2) -> $LOG_FILE"
  else
    echo "Submitting $model_name ..."
    nohup "$RSCRIPT" "$SCRIPT" > "$LOG_FILE" 2>&1 &
    ((running_jobs++)) || true
    if ((running_jobs >= MAX_JOBS)); then
      # Wait for at least one job to finish before starting more
      wait -n || true
      ((running_jobs--)) || true
    fi
  fi
done < <(tail -n +2 "$CONFIG")

# Wait for all remaining background jobs
if [[ $DRY_RUN -eq 0 ]]; then
  wait || true
  echo ""
  echo "All MVMR models finished. Monitor with: tail -f $LOG_DIR/MVMR.${PHASE}.*.log"
fi
