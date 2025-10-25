#!/usr/bin/env bash
set -euo pipefail

## ---------------- Section 0: User-configurable parameters ---------------- ##
BASE_ROOT="path_to_crcount"        # Parent directory containing Cell Ranger outputs
OUT_DIR="path_to_cellbender_out"   # Output directory for CellBender results
THREADS=32                         # CPU threads for CellBender (PyTorch CPU mode)
CHECKPOINT_MIN=120                 # Minutes between periodic checkpoints
CONDA_ENV="cellbender"             # Conda environment name containing CellBender
CONDA_INIT="path_to_conda.sh"      # e.g., "$HOME/miniforge3/etc/profile.d/conda.sh"
EMAIL="your_email@example.com"     # Notification email (requires `mail` on system)

mkdir -p "$OUT_DIR"

## ---------------- Section 1: Global log ---------------------------------- ##
MASTER_LOG="$OUT_DIR/cellbender_batch.log"   # No date suffix (append mode)
{
  echo "=== CellBender batch started $(date) ==="
} | tee -a "$MASTER_LOG"

## ---------------- Section 2: Activate Conda env -------------------------- ##
if [[ -f "$CONDA_INIT" ]]; then
  # shellcheck disable=SC1090
  source "$CONDA_INIT"
  conda activate "$CONDA_ENV"
else
  echo "ERROR: Conda init file not found at: $CONDA_INIT" | tee -a "$MASTER_LOG"
  exit 1
fi
echo "Activated conda env: $(conda info --envs | awk '/\*/{print $1}')" | tee -a "$MASTER_LOG"

## ---------------- Section 3: Iterate samples ----------------------------- ##
# Find all raw_feature_bc_matrix directories under BASE_ROOT
find "$BASE_ROOT" -type d -path "*/outs/raw_feature_bc_matrix" \
| sort | while read -r INPUT; do
  # Sample directory is two levels up from raw_feature_bc_matrix
  SAMPLE_DIR="$(basename "$(dirname "$(dirname "$INPUT")")")"
  SAMPLE="${SAMPLE_DIR%}"
  SAMPLE_OUT="$OUT_DIR/${SAMPLE}"
  OUT_H5="$SAMPLE_OUT/${SAMPLE}.h5"
  LOG_FILE="$SAMPLE_OUT/${SAMPLE}.log"

  mkdir -p "$SAMPLE_OUT"

  echo ">>> [$SAMPLE] START $(date)" | tee -a "$MASTER_LOG" "$LOG_FILE"

  # Run CellBender remove-background
  cellbender remove-background \
    --cpu-threads "$THREADS" \
    --checkpoint-mins "$CHECKPOINT_MIN" \
    --input  "$INPUT" \
    --output "$OUT_H5" \
    >> "$LOG_FILE" 2>&1

  # Move checkpoint archive (if created) into a sample-specific file
  CKPT_SRC="$(dirname "$OUT_H5")/ckpt.tar.gz"
  if [[ -f "$CKPT_SRC" ]]; then
    mv "$CKPT_SRC" "$SAMPLE_OUT/${SAMPLE}_ckpt.tar.gz"
  fi

  echo ">>> [$SAMPLE] DONE  $(date)" | tee -a "$MASTER_LOG" "$LOG_FILE"
done

echo "=== All samples finished $(date) ===" | tee -a "$MASTER_LOG"

## ---------------- Section 4: Email notification -------------------------- ##
if command -v mail >/dev/null 2>&1; then
  {
    echo "CellBender batch has completed."
    echo "Master log: $MASTER_LOG"
  } | mail -s "CellBender batch completed" "$EMAIL"
  echo "Notification email sent to $EMAIL" | tee -a "$MASTER_LOG"
else
  echo "WARNING: 'mail' command not found; no notification sent." | tee -a "$MASTER_LOG"
fi
