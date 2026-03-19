#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
CONFIG_PATH="$PROJECT_ROOT/conf/Config.yaml"
source "$PROJECT_ROOT/script/console_ui.sh"
ui_init
ui_logo

usage() {
    cat <<EOF
Usage:
  bash pipe/run_pipeline.sh [--config conf/Config.yaml]

Options:
  --config PATH   Path to pipeline config YAML
  --help, -h      Show this help message
EOF
}

run_stage() {
    local stage_step="$1"
    local stage_title="$2"
    shift 2
    ui_stage_start "$stage_step" "$stage_title"
    "$@"
    ui_stage_end "$stage_step" "$stage_title"
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            CONFIG_PATH="$2"
            shift 2
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            ui_error "参数错误 | Unknown argument: $1"
            usage >&2
            exit 1
            ;;
    esac
done

if [[ ! -f "$CONFIG_PATH" ]]; then
    ui_error "配置缺失 | Config file not found: $CONFIG_PATH"
    exit 1
fi

ui_section "mtDNA Haplogroup Pipeline" "Build tree index -> annotate samples"
ui_kv "Project root" "$PROJECT_ROOT"
ui_kv "Config" "$CONFIG_PATH"

run_stage "Step 1/4" "environment | 检查运行环境" \
    bash "$PROJECT_ROOT/script/check_env.sh"

ui_stage_start "Step 2/4" "config | 加载配置并解析路径"
source "$PROJECT_ROOT/script/load_config.sh" "$CONFIG_PATH"

require_var() {
    local var_name="$1"
    if [[ -z "${!var_name:-}" ]]; then
        ui_error "配置错误 | Missing required config value: $var_name"
        exit 1
    fi
}

resolve_path() {
    local raw_path="$1"
    if [[ "$raw_path" = /* ]]; then
        printf '%s\n' "$raw_path"
    else
        printf '%s/%s\n' "$PROJECT_ROOT" "$raw_path"
    fi
}

is_truthy() {
    local value="${1,,}"
    [[ "$value" == "1" || "$value" == "true" || "$value" == "yes" ]]
}

require_var CFG_TOOLS_PYTHON
require_var CFG_PATHS_INPUT_FILE
require_var CFG_PATHS_OUTPUT_DIR
require_var CFG_PATHS_TREE_INDEX_FILE
require_var CFG_PATHS_PHYLOTREE_TABLE
require_var CFG_PATHS_INPUT_CORRECTION
require_var CFG_PATHS_LLT_TARGETS
require_var CFG_PATHS_LLT_CORRECTION
require_var CFG_PATHS_YUCHUNLI_TARGETS
require_var CFG_PATHS_YUCHUNLI_CORRECTION
require_var CFG_RUNTIME_OUTPUT_PREFIX
require_var CFG_RUNTIME_REBUILD_INDEX
require_var CFG_RUNTIME_LOG_LEVEL

PYTHON_BIN="$CFG_TOOLS_PYTHON"
INPUT_FILE=$(resolve_path "$CFG_PATHS_INPUT_FILE")
OUTPUT_DIR=$(resolve_path "$CFG_PATHS_OUTPUT_DIR")
TREE_INDEX_PATH=$(resolve_path "$CFG_PATHS_TREE_INDEX_FILE")
PHYLOTREE_TABLE=$(resolve_path "$CFG_PATHS_PHYLOTREE_TABLE")
INPUT_CORRECTION=$(resolve_path "$CFG_PATHS_INPUT_CORRECTION")
LLT_TARGETS=$(resolve_path "$CFG_PATHS_LLT_TARGETS")
LLT_CORRECTION=$(resolve_path "$CFG_PATHS_LLT_CORRECTION")
YUCHUNLI_TARGETS=$(resolve_path "$CFG_PATHS_YUCHUNLI_TARGETS")
YUCHUNLI_CORRECTION=$(resolve_path "$CFG_PATHS_YUCHUNLI_CORRECTION")
OUTPUT_PREFIX="$CFG_RUNTIME_OUTPUT_PREFIX"
REBUILD_INDEX="$CFG_RUNTIME_REBUILD_INDEX"
LOG_LEVEL="$CFG_RUNTIME_LOG_LEVEL"

mkdir -p "$OUTPUT_DIR"
mkdir -p "$(dirname "$TREE_INDEX_PATH")"
ANNOTATION_PATH="$OUTPUT_DIR/${OUTPUT_PREFIX}.tsv"

ui_kv "Input" "$INPUT_FILE"
ui_kv "Tree index" "$TREE_INDEX_PATH"
ui_kv "Output" "$ANNOTATION_PATH"
ui_kv "Rebuild index" "$REBUILD_INDEX"
ui_kv "Log level" "$LOG_LEVEL"
ui_stage_end "Step 2/4" "config | 加载配置并解析路径"

if is_truthy "$REBUILD_INDEX" || [[ ! -f "$TREE_INDEX_PATH" ]]; then
    run_stage "Step 3/4" "tree index | 构建或刷新索引" \
        "$PYTHON_BIN" "$PROJECT_ROOT/python/build_tree_index.py" \
        --phylotree-table "$PHYLOTREE_TABLE" \
        --input-correction "$INPUT_CORRECTION" \
        --llt-targets "$LLT_TARGETS" \
        --llt-correction "$LLT_CORRECTION" \
        --yuchunli-targets "$YUCHUNLI_TARGETS" \
        --yuchunli-correction "$YUCHUNLI_CORRECTION" \
        --output "$TREE_INDEX_PATH" \
        --log-level "$LOG_LEVEL"
else
    ui_stage_start "Step 3/4" "tree index | 复用已有索引"
    ui_info "Reusing existing tree index: $TREE_INDEX_PATH"
    ui_stage_end "Step 3/4" "tree index | 复用已有索引"
fi

run_stage "Step 4/4" "annotate | 样本单倍群注释" \
    "$PYTHON_BIN" "$PROJECT_ROOT/python/annotate_haplogroups.py" \
    --input "$INPUT_FILE" \
    --tree-index "$TREE_INDEX_PATH" \
    --output "$ANNOTATION_PATH" \
    --log-level "$LOG_LEVEL"

ui_section "Workflow completed" "mtDNA haplogroup annotation finished successfully."
ui_summary_line "Tree index" "$TREE_INDEX_PATH"
ui_summary_line "Result" "$ANNOTATION_PATH"
