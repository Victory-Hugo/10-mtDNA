#!/usr/bin/env bash
set -euo pipefail

CONFIG_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/13-phylotree-update/1-excel→json/conf/Config.yaml"
source "/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/13-phylotree-update/1-excel→json/script/load_config.sh" "$CONFIG_FILE"

INPUT_EXCEL="$(resolve_path "$PATHS_INPUT_EXCEL")"
TEMPLATE_JSON="$(resolve_path "$PATHS_TEMPLATE_JSON")"
OUTPUT_JSON="$(resolve_path "$PATHS_OUTPUT_JSON")"
OUTPUT_TSV="$(resolve_path "$PATHS_OUTPUT_TSV")"
PYTHON_MODULE="$(resolve_path "$TOOLS_PYTHON_MODULE")"
LOG_DIR="$(resolve_path "$PATHS_LOG_DIR")"
OUTPUT_DIR="$(resolve_path "$PATHS_OUTPUT_DIR")"
TEMP_DIR="$(resolve_path "$PATHS_TEMP_DIR")"

mkdir -p "$LOG_DIR" "$OUTPUT_DIR" "$TEMP_DIR"

if [[ "$RUNTIME_OVERWRITE" != "true" && -e "$OUTPUT_JSON" ]]; then
  echo "输出文件已存在且 overwrite=false: $OUTPUT_JSON" >&2
  exit 1
fi

"$TOOLS_CONDA_BIN" run -n "$TOOLS_CONDA_ENV_NAME" python "$PYTHON_MODULE" \
  --input-excel "$INPUT_EXCEL" \
  --template-json "$TEMPLATE_JSON" \
  --output-json "$OUTPUT_JSON" \
  --output-tsv "$OUTPUT_TSV" \
  --sheet-name "$RUNTIME_PHYLOTREE_SHEET" \
  --correction-sheet "$RUNTIME_CORRECTION_SHEET" \
  --indent "$RUNTIME_INDENT" \
  2>&1 | tee "$LOG_DIR/1-excel_to_json.log"
