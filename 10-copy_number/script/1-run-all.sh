#!/usr/bin/bash
set -euo pipefail

BASE_DIR="$(cd "$(dirname "$0")/.." && pwd)"
CONFIG="$BASE_DIR/conf/config.yaml"
JOBS=""
FORCE=0

read_config() {
  local config_path="$1"
  local key="$2"
  local default="${3-__MISSING__}"
  local value

  value=$(awk -v key="$key" '
    /^[[:space:]]*#/ {next}
    /^[[:space:]]*$/ {next}
    /^[^[:space:]][^:]*:[[:space:]]*/ {
      k=$0
      sub(/:.*/, "", k)
      gsub(/^[[:space:]]+|[[:space:]]+$/, "", k)
      if (k == key) {
        v=$0
        sub(/^[^:]*:[[:space:]]*/, "", v)
        print v
        exit 0
      }
    }
  ' "$config_path")

  if [[ -z "$value" ]]; then
    if [[ "$default" == "__MISSING__" ]]; then
      echo "Missing key in config: $key" >&2
      return 1
    fi
    echo "$default"
    return 0
  fi

  if [[ "$value" =~ ^[\{\[] ]]; then
    echo "Key is not a scalar value: $key" >&2
    return 1
  fi

  value="${value%%[[:space:]]#*}"
  value="${value#"${value%%[![:space:]]*}"}"
  value="${value%"${value##*[![:space:]]}"}"
  if [[ "$value" =~ ^\".*\"$ ]]; then
    value="${value:1:${#value}-2}"
  elif [[ "$value" =~ ^\'.*\'$ ]]; then
    value="${value:1:${#value}-2}"
  fi

  echo "$value"
}

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

LOG_DIR=$(read_config "$CONFIG" log_dir)
SUCCESS_LOG=$(read_config "$CONFIG" success_log)
LIST_PATH=$(read_config "$CONFIG" list_path)
INPUT_TYPE=$(read_config "$CONFIG" input_type "BAM")
MT_COPY_PY=$(read_config "$CONFIG" mt_copy_number_py)
MERGE_PY=$(read_config "$CONFIG" merge_results_py)
MERGE_NAME=$(read_config "$CONFIG" merge_output_name)
OUTPUT_DIR=$(read_config "$CONFIG" output_dir)
JOBS_DEFAULT=$(read_config "$CONFIG" jobs_default)

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
    if [[ -s "$SUCCESS_LOG" ]]; then
      awk 'FNR==NR {done[$1]=1; next} !($1 in done)' "$SUCCESS_LOG" "$LIST_PATH" | \
        parallel --colsep '\t' --jobs "$JOBS" --halt now,fail=1 \
          python "$MT_COPY_PY" \
            --config "$CONFIG" --sample-id {1} --input {2} --input-type "$INPUT_TYPE" $force_arg --success-log "$SUCCESS_LOG"
    else
      parallel --colsep '\t' --jobs "$JOBS" --halt now,fail=1 \
        python "$MT_COPY_PY" \
          --config "$CONFIG" --sample-id {1} --input {2} --input-type "$INPUT_TYPE" $force_arg --success-log "$SUCCESS_LOG" \
        :::: "$LIST_PATH"
    fi
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
    if [[ -s "$SUCCESS_LOG" ]]; then
      awk 'FNR==NR {done[$1]=1; next} !($1 in done)' "$SUCCESS_LOG" "$LIST_PATH" | \
        xargs -P "$JOBS" -n 2 sh -c \
          'python "$MT_COPY_PY" --config "$CONFIG" --sample-id "$1" --input "$2" --input-type "$INPUT_TYPE" $FORCE_ARG --success-log "$SUCCESS_LOG"' \
          _
    else
      cut -f1,2 "$LIST_PATH" | xargs -P "$JOBS" -n 2 sh -c \
        'python "$MT_COPY_PY" --config "$CONFIG" --sample-id "$1" --input "$2" --input-type "$INPUT_TYPE" $FORCE_ARG --success-log "$SUCCESS_LOG"' \
        _
    fi
  fi
fi

python "$MERGE_PY" \
  --output-dir "$OUTPUT_DIR" \
  --output "$OUTPUT_DIR/$MERGE_NAME"
