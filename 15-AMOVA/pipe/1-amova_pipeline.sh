#!/usr/bin/env bash
# =============================================================================
# AMOVA Pipeline — 总控脚本
# 用法：bash pipe/1-amova_pipeline.sh [config_path]
# 默认 config_path = <项目根目录>/conf/Config.yaml
#
# 断点续跑机制
# ─────────────────────────────────────────────────────────────────────────────
# Step 1（VCF 解析）：计算输入指纹（VCF 大小+时间戳、分组文件 MD5、参数串）。
#   若指纹与上次相同且所有输出文件完整，则跳过，耗时 < 1s。
#
# Step 2（AMOVA 计算）：每个场景独立维护缓存键（场景定义 + 置换参数）。
#   overwrite=false 时：指纹匹配则跳过该场景；
#   overwrite=true  时：始终重算，但同时刷新缓存键。
#
# 要强制全量重跑，设 overwrite: true 并删除各 .amova_cache 文件，
# 或直接 rm -rf output/ temp/。
# =============================================================================
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "$SCRIPT_DIR/.." && pwd)
CONFIG_PATH="${1:-$PROJECT_ROOT/conf/Config.yaml}"

# ── UI 库（在日志重定向前 source + init，保证终端检测正确）─────────────────
source "$PROJECT_ROOT/script/amova_ui.sh"
amova_ui_init

# ── 配置文件检查 ──────────────────────────────────────────────────────────────
if [[ ! -f "$CONFIG_PATH" ]]; then
    ui_error "配置文件不存在: $CONFIG_PATH"
    exit 1
fi

ui_banner

# ── 加载配置 ──────────────────────────────────────────────────────────────────
ui_info "加载配置: $CONFIG_PATH"
source "$PROJECT_ROOT/script/load_config.sh" "$CONFIG_PATH"

for var in CFG_PROJECT_BASE_DIR CFG_PROJECT_INPUT_VCF CFG_PROJECT_INPUT_SCENARIOS \
           CFG_PROJECT_CONTIG_NAME CFG_TOOLS_PYTHON_BIN \
           CFG_PATHS_PYTHON_DIR CFG_PATHS_LOG_DIR \
           CFG_PATHS_OUTPUT_STEP1_TABLE CFG_PATHS_TEMP_STEP1 \
           CFG_PATHS_OUTPUT_STEP2_TABLE CFG_PATHS_OUTPUT_STEP2_REPORT \
           CFG_RUNTIME_MIN_POP_SIZE CFG_RUNTIME_PERMUTATION_N \
           CFG_RUNTIME_RANDOM_SEED CFG_RUNTIME_OVERWRITE; do
    if [[ -z "${!var:-}" ]]; then
        ui_error "配置缺少必需变量: $var"
        exit 1
    fi
done

# ── 路径解析 ──────────────────────────────────────────────────────────────────
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

REF_GROUP_FILE=$(resolve_path "input/1-group.tsv")

# ── 配置摘要 ──────────────────────────────────────────────────────────────────
ui_kv "Input VCF"      "$(basename "$INPUT_VCF")"
ui_kv "Scenarios"      "$(basename "$SCENARIOS_YAML")"
ui_kv "Contig"         "$CONTIG_NAME"
ui_kv "Min pop size"   "$MIN_POP_SIZE"
ui_kv "Permutations"   "$PERM_N"
ui_kv "Random seed"    "$RAND_SEED"
ui_kv "Parallel jobs"  "$N_JOBS"
ui_kv "Overwrite"      "$OVERWRITE"

# ── 建立输出目录 ──────────────────────────────────────────────────────────────
mkdir -p "$LOG_DIR" "$STEP1_TABLE" "$STEP1_TEMP" "$STEP2_TABLE" "$STEP2_REPORT"

# ── 日志重定向 ────────────────────────────────────────────────────────────────
LOG_FILE="$LOG_DIR/pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "$LOG_FILE") 2>&1
ui_info "日志: $LOG_FILE"

# ── 环境检查 ──────────────────────────────────────────────────────────────────
bash "$PROJECT_ROOT/script/check_env.sh" "$CONFIG_PATH"

T_PIPELINE=$SECONDS

# ══════════════════════════════════════════════════════════════════════════════
# 缓存指纹工具函数
# ══════════════════════════════════════════════════════════════════════════════

# 返回文件"大小_修改时间"，不读取文件内容（对大 VCF 极快）
_file_stat() {
    stat -c "%s_%Y" "$1" 2>/dev/null \
        || stat -f "%z_%m" "$1" 2>/dev/null \
        || printf 'nohash'
}

# 对字符串计算 MD5（用于拼接多个字段形成唯一指纹）
_md5_str() {
    printf '%s' "$1" | md5sum | cut -d' ' -f1
}

# ══════════════════════════════════════════════════════════════════════════════
# Step 1  VCF 解析 → 频率矩阵 + 基因型矩阵
# ══════════════════════════════════════════════════════════════════════════════

# 指纹组成：VCF 大小+时间戳 | 分组文件 MD5 | 染色体名 | 最小种群阈值
# 任意一项改变都会触发重新解析
_vcf_stat=$(_file_stat "$INPUT_VCF")
_grp_md5=$(md5sum "$REF_GROUP_FILE" | cut -d' ' -f1)
STEP1_KEY=$(_md5_str "${_vcf_stat}|${_grp_md5}|${CONTIG_NAME}|${MIN_POP_SIZE}")

STEP1_CACHE_FILE="$STEP1_TABLE/.amova_cache"

_step1_outputs_exist() {
    [[ -f "$STEP1_TABLE/allele_freq_matrix.npz" ]] && \
    [[ -f "$STEP1_TABLE/pop_names.txt"          ]] && \
    [[ -f "$STEP1_TABLE/sample_group_map.tsv"   ]] && \
    [[ -f "$STEP1_TABLE/vcf_parse_summary.tsv"  ]]
}

_step1_cache_valid() {
    [[ "$OVERWRITE" != "true" ]] && \
    [[ -f "$STEP1_CACHE_FILE" ]] && \
    [[ "$(cat "$STEP1_CACHE_FILE" 2>/dev/null)" == "$STEP1_KEY" ]] && \
    _step1_outputs_exist
}

if _step1_cache_valid; then
    ui_cache_hit 1 "VCF 解析" "输入未变更，跳过重新解析"
else
    if [[ "$OVERWRITE" == "true" ]]; then
        ui_cache_miss "overwrite=true，强制重新解析"
    elif ! _step1_outputs_exist; then
        ui_cache_miss "输出文件缺失，重新解析"
    else
        ui_cache_miss "输入已变更，缓存失效，重新解析"
    fi

    ui_step_begin 1 "VCF 解析"
    T1=$SECONDS

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

    # 写入缓存键（仅在解析成功后）
    printf '%s' "$STEP1_KEY" > "$STEP1_CACHE_FILE"
    ui_step_done 1 "VCF 解析" $(( SECONDS - T1 ))
fi

# ══════════════════════════════════════════════════════════════════════════════
# Step 2  AMOVA 计算（所有场景；逐场景缓存由 Python 端处理）
# ══════════════════════════════════════════════════════════════════════════════
ui_step_begin 2 "AMOVA 计算"
T2=$SECONDS

"$PYTHON_BIN" "$PYTHON_DIR/1-3-run_amova.py" \
    --freq-matrix    "$STEP1_TABLE/allele_freq_matrix.npz" \
    --pop-names      "$STEP1_TABLE/pop_names.txt" \
    --scenarios      "$SCENARIOS_YAML" \
    --base-dir       "$BASE_DIR" \
    --output-dir     "$STEP2_TABLE" \
    --permutation-n  "$PERM_N" \
    --random-seed    "$RAND_SEED" \
    --overwrite      "$([ "$OVERWRITE" = "true" ] && echo 1 || echo 0)" \
    --n-jobs         "$N_JOBS" \
    --step1-key      "$STEP1_KEY"

ui_step_done 2 "AMOVA 计算" $(( SECONDS - T2 ))

# ══════════════════════════════════════════════════════════════════════════════
# Step 3  汇总结果 → Markdown 报告
# ══════════════════════════════════════════════════════════════════════════════
ui_step_begin 3 "生成报告"
T3=$SECONDS

"$PYTHON_BIN" "$PYTHON_DIR/1-4-collect_results.py" \
    --results-dir   "$STEP2_TABLE" \
    --scenarios     "$SCENARIOS_YAML" \
    --output-md     "$STEP2_REPORT/amova_summary.md"

ui_step_done 3 "生成报告" $(( SECONDS - T3 ))

# ══════════════════════════════════════════════════════════════════════════════
# 完成摘要
# ══════════════════════════════════════════════════════════════════════════════
ui_final_summary
ui_stat "Report"        "$STEP2_REPORT/amova_summary.md"
ui_stat "Log"           "$LOG_FILE"
ui_stat "Total elapsed" "$(( SECONDS - T_PIPELINE ))s"
printf '\n' >&2
