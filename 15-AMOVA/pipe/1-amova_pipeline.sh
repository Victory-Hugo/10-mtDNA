#!/usr/bin/env bash
# =============================================================================
# AMOVA Pipeline — 总控脚本
# 用法：bash pipe/1-amova_pipeline.sh [config_path]
# 默认 config_path = <项目根目录>/conf/Config.yaml
# =============================================================================
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
CONFIG_PATH="${1:-$PROJECT_ROOT/conf/Config.yaml}"

if [[ ! -f "$CONFIG_PATH" ]]; then
    echo "[ERROR] 配置文件不存在: $CONFIG_PATH" >&2
    exit 1
fi

# ── 1. 加载配置 ──────────────────────────────────────────────────────────────
source "$PROJECT_ROOT/script/load_config.sh" "$CONFIG_PATH"

# 必需变量检查
for var in CFG_PROJECT_BASE_DIR CFG_PROJECT_INPUT_VCF CFG_PROJECT_INPUT_SCENARIOS \
           CFG_PROJECT_CONTIG_NAME CFG_TOOLS_PYTHON_BIN \
           CFG_PATHS_PYTHON_DIR CFG_PATHS_LOG_DIR \
           CFG_PATHS_OUTPUT_STEP1_TABLE CFG_PATHS_TEMP_STEP1 \
           CFG_PATHS_OUTPUT_STEP2_TABLE CFG_PATHS_OUTPUT_STEP2_REPORT \
           CFG_RUNTIME_MIN_POP_SIZE CFG_RUNTIME_PERMUTATION_N \
           CFG_RUNTIME_RANDOM_SEED CFG_RUNTIME_OVERWRITE; do
    if [[ -z "${!var:-}" ]]; then
        echo "[ERROR] 配置缺少变量: $var" >&2
        exit 1
    fi
done

# ── 2. 路径解析（相对路径相对于 PROJECT_ROOT）────────────────────────────────
resolve_path() {
    local raw="$1"
    if [[ "$raw" = /* ]]; then printf '%s\n' "$raw"
    else printf '%s\n' "$PROJECT_ROOT/$raw"; fi
}

BASE_DIR=$(resolve_path "$CFG_PROJECT_BASE_DIR")
INPUT_VCF=$(resolve_path "$CFG_PROJECT_INPUT_VCF")
SCENARIOS_YAML=$(resolve_path "$CFG_PROJECT_INPUT_SCENARIOS")
CONTIG_NAME="$CFG_PROJECT_CONTIG_NAME"
PYTHON_BIN=$(resolve_path "$CFG_TOOLS_PYTHON_BIN")
PYTHON_DIR=$(resolve_path "$CFG_PATHS_PYTHON_DIR")

LOG_DIR=$(resolve_path "$CFG_PATHS_LOG_DIR")
STEP1_TABLE=$(resolve_path "$CFG_PATHS_OUTPUT_STEP1_TABLE")
STEP1_TEMP=$(resolve_path "$CFG_PATHS_TEMP_STEP1")
STEP2_TABLE=$(resolve_path "$CFG_PATHS_OUTPUT_STEP2_TABLE")
STEP2_REPORT=$(resolve_path "$CFG_PATHS_OUTPUT_STEP2_REPORT")

MIN_POP_SIZE="$CFG_RUNTIME_MIN_POP_SIZE"
PERM_N="$CFG_RUNTIME_PERMUTATION_N"
RAND_SEED="$CFG_RUNTIME_RANDOM_SEED"
OVERWRITE="$CFG_RUNTIME_OVERWRITE"
N_JOBS="${CFG_RUNTIME_N_JOBS:-1}"

# ── 3. 建立输出目录 ───────────────────────────────────────────────────────────
mkdir -p "$LOG_DIR" "$STEP1_TABLE" "$STEP1_TEMP" "$STEP2_TABLE" "$STEP2_REPORT"

# ── 4. 日志重定向 ─────────────────────────────────────────────────────────────
LOG_FILE="$LOG_DIR/pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Pipeline 启动"
echo "  项目根目录: $PROJECT_ROOT"
echo "  输入 VCF:   $INPUT_VCF"
echo "  场景配置:   $SCENARIOS_YAML"

# ── 5. 环境检查 ───────────────────────────────────────────────────────────────
bash "$PROJECT_ROOT/script/check_env.sh" "$CONFIG_PATH"

# ══════════════════════════════════════════════════════════════════════════════
# 步骤 1：解析 VCF → 频率矩阵 + 基因型矩阵
# ══════════════════════════════════════════════════════════════════════════════
echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] ── 步骤1：VCF 解析 ──"

# 使用 1-group.tsv 作为参考分组文件（含两级分组信息，最完整）
REF_GROUP_FILE=$(resolve_path "input/1-group.tsv")

"$PYTHON_BIN" "$PYTHON_DIR/1-1-parse_vcf.py" \
    --vcf                 "$INPUT_VCF" \
    --group-file          "$REF_GROUP_FILE" \
    --id-col              "ID" \
    --pop-col             "Group_small" \
    --contig-name         "$CONTIG_NAME" \
    --min-pop-size        "$MIN_POP_SIZE" \
    --output-freq-matrix  "$STEP1_TABLE/allele_freq_matrix.npz" \
    --output-pop-names    "$STEP1_TABLE/pop_names.txt" \
    --output-sample-map   "$STEP1_TABLE/sample_group_map.tsv" \
    --output-summary      "$STEP1_TABLE/vcf_parse_summary.tsv"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] 步骤1 完成"

# ══════════════════════════════════════════════════════════════════════════════
# 步骤 2：为所有场景运行 AMOVA
# ══════════════════════════════════════════════════════════════════════════════
echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] ── 步骤2：AMOVA 计算 ──"

"$PYTHON_BIN" "$PYTHON_DIR/1-3-run_amova.py" \
    --freq-matrix    "$STEP1_TABLE/allele_freq_matrix.npz" \
    --pop-names      "$STEP1_TABLE/pop_names.txt" \
    --scenarios      "$SCENARIOS_YAML" \
    --base-dir       "$BASE_DIR" \
    --output-dir     "$STEP2_TABLE" \
    --permutation-n  "$PERM_N" \
    --random-seed    "$RAND_SEED" \
    --overwrite      "$([ "$OVERWRITE" = "true" ] && echo 1 || echo 0)" \
    --n-jobs         "$N_JOBS"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] 步骤2 完成"

# ══════════════════════════════════════════════════════════════════════════════
# 步骤 3：汇总结果 → Markdown 报告
# ══════════════════════════════════════════════════════════════════════════════
echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] ── 步骤3：生成报告 ──"

"$PYTHON_BIN" "$PYTHON_DIR/1-4-collect_results.py" \
    --results-dir   "$STEP2_TABLE" \
    --scenarios     "$SCENARIOS_YAML" \
    --output-md     "$STEP2_REPORT/amova_summary.md"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] 步骤3 完成"

# ── 完成 ──────────────────────────────────────────────────────────────────────
echo ""
echo "[$(date '+%Y-%m-%d %H:%M:%S')] Pipeline 完成"
echo "  结果报告: $STEP2_REPORT/amova_summary.md"
echo "  日志文件: $LOG_FILE"
