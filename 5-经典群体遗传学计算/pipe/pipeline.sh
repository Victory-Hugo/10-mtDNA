#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
#! 修改CONFIG_PATH为实际配置文件路径
CONFIG_PATH="${1:-$ROOT_DIR/conf/Config.bootstrap_demo.json}" 

if [[ ! -f "$CONFIG_PATH" ]]; then
  echo "Config file not found: $CONFIG_PATH" >&2
  exit 1
fi

eval "$(python3 "$ROOT_DIR/python/00_emit_config_env.py" --config "$CONFIG_PATH")"

mkdir -p "$OUTPUT_DIR" "$TMP_DIR" "$META_DIR" "$LOG_DIR" \
  "$OUTPUT_DIR/0-qc" "$OUTPUT_DIR/1-within" "$OUTPUT_DIR/2-between" "$OUTPUT_DIR/3-windows"

LOG_FILE="$LOG_DIR/pipeline.log"
exec > >(tee -a "$LOG_FILE") 2>&1

echo "[$(date '+%F %T')] Pipeline started"
echo "Project root: $PROJECT_ROOT"
echo "Config: $CONFIG_PATH"

"$PYTHON_BIN" "$ROOT_DIR/python/01_prepare_inputs.py" \
  --input-vcf "$INPUT_VCF" \
  --sample-table "$SAMPLE_TABLE" \
  --sample-id-column "$SAMPLE_ID_COLUMN" \
  --group-column "$GROUP_COLUMN" \
  --min-group-size "$MIN_GROUP_SIZE" \
  --matched-samples-output "$MATCHED_SAMPLES_TSV" \
  --group-counts-output "$GROUP_COUNTS_TSV" \
  --summary-output "$INPUT_SUMMARY_TSV"

"$PYTHON_BIN" "$ROOT_DIR/python/02_vcf_to_zarr.py" \
  --input-vcf "$INPUT_VCF" \
  --matched-samples "$MATCHED_SAMPLES_TSV" \
  --sample-id-column "$SAMPLE_ID_COLUMN" \
  --zarr-path "$ZARR_PATH" \
  --overwrite "$OVERWRITE_ZARR" \
  --chunk-length "$ZARR_CHUNK_LENGTH" \
  --chunk-width "$ZARR_CHUNK_WIDTH" \
  --alt-number "$ALT_NUMBER"

"$PYTHON_BIN" "$ROOT_DIR/python/03_compute_diversity.py" \
  --zarr-path "$ZARR_PATH" \
  --matched-samples "$MATCHED_SAMPLES_TSV" \
  --sample-id-column "$SAMPLE_ID_COLUMN" \
  --group-column "$GROUP_COLUMN" \
  --contig-name "$CONTIG_NAME" \
  --sequence-length "$SEQUENCE_LENGTH" \
  --within-metrics "$WITHIN_METRICS" \
  --between-metrics "$BETWEEN_METRICS" \
  --window-enable "$WINDOW_ENABLE" \
  --window-size "$WINDOW_SIZE" \
  --window-step "$WINDOW_STEP" \
  --within-output "$WITHIN_SUMMARY_TSV" \
  --between-output "$BETWEEN_SUMMARY_TSV" \
  --window-output "$WINDOW_METRICS_TSV"

if [[ "$RESAMPLING_ENABLE" == "1" ]]; then
  mkdir -p "$BOOTSTRAP_DIR"
  "$PYTHON_BIN" "$ROOT_DIR/python/04_bootstrap_diversity.py" \
    --zarr-path "$ZARR_PATH" \
    --matched-samples "$MATCHED_SAMPLES_TSV" \
    --sample-id-column "$SAMPLE_ID_COLUMN" \
    --group-column "$GROUP_COLUMN" \
    --contig-name "$CONTIG_NAME" \
    --sequence-length "$SEQUENCE_LENGTH" \
    --within-metrics "$WITHIN_METRICS" \
    --between-metrics "$BETWEEN_METRICS" \
    --bootstrap-replicates "$BOOTSTRAP_REPLICATES" \
    --random-seed "$BOOTSTRAP_RANDOM_SEED" \
    --sample-size-strategy "$RESAMPLE_STRATEGY" \
    --within-reference "$WITHIN_SUMMARY_TSV" \
    --between-reference "$BETWEEN_SUMMARY_TSV" \
    --within-bootstrap-output "$WITHIN_BOOTSTRAP_SUMMARY_TSV" \
    --between-bootstrap-output "$BETWEEN_BOOTSTRAP_SUMMARY_TSV" \
    --run-summary-output "$BOOTSTRAP_RUN_SUMMARY_TSV" \
    --write-replicate-tables "$WRITE_BOOTSTRAP_REPLICATE_TABLES" \
    --within-replicates-output "$WITHIN_BOOTSTRAP_REPLICATES_TSV" \
    --between-replicates-output "$BETWEEN_BOOTSTRAP_REPLICATES_TSV" \
    --n-threads "$N_THREADS"
else
  echo "Bootstrap disabled by config"
fi

echo "[$(date '+%F %T')] Pipeline finished"
