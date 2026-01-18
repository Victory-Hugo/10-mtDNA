#!/usr/bin/bash
set -euo pipefail

BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
CONFIG="$BASE_DIR/conf/config.yaml"
JOBS=""
FORCE=0

while [[ $# -gt 0 ]]; do
  case "$1" in
    --config)
      CONFIG="$2"
      shift 2
      ;;
    --jobs)
      JOBS="$2"
      shift 2
      ;;
    --force)
      FORCE=1
      shift
      ;;
    *)
      echo "Unknown argument: $1" >&2
      exit 1
      ;;
  esac
done

CONFIG_CLI="$BASE_DIR/python/config_cli.py"
LOG_DIR=$(python "$CONFIG_CLI" --config "$CONFIG" --key log_dir)
SUCCESS_LOG=$(python "$CONFIG_CLI" --config "$CONFIG" --key success_log)
LIST_PATH=$(python "$CONFIG_CLI" --config "$CONFIG" --key list_path)
INPUT_TYPE=$(python "$CONFIG_CLI" --config "$CONFIG" --key input_type --default "BAM")
MT_COPY_PY=$(python "$CONFIG_CLI" --config "$CONFIG" --key mt_copy_number_py)
MERGE_PY=$(python "$CONFIG_CLI" --config "$CONFIG" --key merge_results_py)
MERGE_NAME=$(python "$CONFIG_CLI" --config "$CONFIG" --key merge_output_name)
OUTPUT_DIR=$(python "$CONFIG_CLI" --config "$CONFIG" --key output_dir)
JOBS_DEFAULT=$(python "$CONFIG_CLI" --config "$CONFIG" --key jobs_default)

if [[ -z "$JOBS" ]]; then
  JOBS="$JOBS_DEFAULT"
fi

mkdir -p "$LOG_DIR"
if [[ -n "$SUCCESS_LOG" ]]; then
  mkdir -p "$(dirname "$SUCCESS_LOG")"
  touch "$SUCCESS_LOG"
fi
LOG_FILE="$LOG_DIR/1-run-all.sh.log"

if [[ -f "$LOG_FILE" ]] && grep -q "STATUS\tSUCCESS" "$LOG_FILE" && [[ "$FORCE" -eq 0 ]]; then
  echo -e "STATUS\tSKIP\tprevious success" >> "$LOG_FILE"
  exit 0
fi

exec > >(tee -a "$LOG_FILE") 2>&1

status=0
trap 'status=$?; if [[ $status -ne 0 ]]; then echo -e "STATUS\tFAILED\tcode=$status"; else echo -e "STATUS\tSUCCESS"; fi' EXIT

echo "CONFIG\t$CONFIG"
echo "LIST\t$LIST_PATH"
echo "INPUT_TYPE\t$INPUT_TYPE"
echo "OUTPUT_DIR\t$OUTPUT_DIR"
echo "SUCCESS_LOG\t$SUCCESS_LOG"

echo "JOB_MODE\tparallel"

force_arg=""
if [[ "$FORCE" -eq 1 ]]; then
  force_arg="--force"
fi

if command -v parallel >/dev/null 2>&1; then
  if [[ "$FORCE" -eq 1 ]]; then
    parallel --colsep '\t' --jobs "$JOBS" --halt now,fail=1 \
      python "$MT_COPY_PY" \
        --config "$CONFIG" --sample-id {1} --input {2} --input-type "$INPUT_TYPE" $force_arg --success-log "$SUCCESS_LOG" \
      :::: "$LIST_PATH"
  else
    awk 'FNR==NR {done[$1]=1; next} !($1 in done)' "$SUCCESS_LOG" "$LIST_PATH" | \
      parallel --colsep '\t' --jobs "$JOBS" --halt now,fail=1 \
        python "$MT_COPY_PY" \
          --config "$CONFIG" --sample-id {1} --input {2} --input-type "$INPUT_TYPE" $force_arg --success-log "$SUCCESS_LOG"
  fi
else
  echo "GNU parallel not found, falling back to xargs -P" >&2
  export CONFIG
  export FORCE_ARG="$force_arg"
  export SUCCESS_LOG
  export MT_COPY_PY
  export INPUT_TYPE
  if [[ "$FORCE" -eq 1 ]]; then
    cut -f1,2 "$LIST_PATH" | xargs -P "$JOBS" -n 2 sh -c \
      'python "$MT_COPY_PY" --config "$CONFIG" --sample-id "$1" --input "$2" --input-type "$INPUT_TYPE" $FORCE_ARG --success-log "$SUCCESS_LOG"' \
      _
  else
    awk 'FNR==NR {done[$1]=1; next} !($1 in done)' "$SUCCESS_LOG" "$LIST_PATH" | \
      xargs -P "$JOBS" -n 2 sh -c \
        'python "$MT_COPY_PY" --config "$CONFIG" --sample-id "$1" --input "$2" --input-type "$INPUT_TYPE" $FORCE_ARG --success-log "$SUCCESS_LOG"' \
        _
  fi
fi

python "$MERGE_PY" \
  --output-dir "$OUTPUT_DIR" \
  --output "$OUTPUT_DIR/$MERGE_NAME"
