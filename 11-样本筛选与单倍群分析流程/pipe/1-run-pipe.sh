#!/bin/bash
# -*- coding: utf-8 -*-
# 样本筛选与单倍群分析流程主管道
# 配置统一从 conf/Config.yaml 读取

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
CONF_DIR="${PROJECT_DIR}/conf"
CONFIG_FILE="${CONF_DIR}/Config.yaml"
LOAD_CONFIG_SH="${PROJECT_DIR}/script/load_config.sh"
UI_SCRIPT="${PROJECT_DIR}/script/console_ui.sh"

if [ ! -f "$CONFIG_FILE" ]; then
    echo "[ERROR] 未找到配置文件: $CONFIG_FILE" >&2
    exit 1
fi
if [ ! -f "$LOAD_CONFIG_SH" ]; then
    echo "[ERROR] 未找到配置加载脚本: $LOAD_CONFIG_SH" >&2
    exit 1
fi
if [ ! -f "$UI_SCRIPT" ]; then
    echo "[ERROR] 未找到控制台 UI 脚本: $UI_SCRIPT" >&2
    exit 1
fi

# shellcheck source=/dev/null
source "$LOAD_CONFIG_SH" "$CONFIG_FILE"
# shellcheck source=/dev/null
source "$UI_SCRIPT"

split_csv_to_array() {
    local raw="$1"
    local -n result_ref=$2
    local item
    result_ref=()
    IFS=',' read -r -a result_ref <<< "$raw"
    for item in "${!result_ref[@]}"; do
        result_ref[$item]="$(echo "${result_ref[$item]}" | sed 's/^[[:space:]]*//; s/[[:space:]]*$//')"
    done
}

normalize_path() {
    local raw_path="$1"
    if [[ "$raw_path" = /* ]]; then
        printf '%s\n' "$raw_path"
    else
        printf '%s/%s\n' "$BASE_DIR" "$raw_path"
    fi
}

# ============================= 配置映射 ============================= #
BASE_DIR="${PROJECT_BASE_DIR}"
INPUT_EXCEL="${PROJECT_INPUT_EXCEL}"
PYTHON3="${TOOLS_PYTHON_BIN}"
RSCRIPT="${TOOLS_RSCRIPT_BIN}"
CONDA_BIN="${TOOLS_CONDA_BIN}"
PYTHON_ENV_PREFIX="${TOOLS_PYTHON_ENV_PREFIX}"

PYTHON_DIR="$(normalize_path "$PATHS_PYTHON_DIR")"
SRC_DIR="$(normalize_path "$PATHS_SRC_DIR")"
HELPER_DIR="$(normalize_path "$PATHS_HELPER_DIR")"
DATA_DIR="$(normalize_path "$PATHS_DATA_DIR")"
TMP_DIR="$(normalize_path "$PATHS_TMP_DIR")"
LOG_DIR="$(normalize_path "$PATHS_LOG_DIR")"
OUTPUT_DIR="$(normalize_path "$PATHS_OUTPUT_DIR")"
PCA_OUTPUT_DIR_GLOBAL="$(normalize_path "$PATHS_PCA_OUTPUT_DIR_GLOBAL")"
PCA_OUTPUT_DIR_CHINA="$(normalize_path "$PATHS_PCA_OUTPUT_DIR_CHINA")"
PCA_VIS_OUTPUT_DIR="$(normalize_path "$PATHS_PCA_VISUALIZATION_DIR")"
HEATMAP_OUTPUT_DIR="$(normalize_path "$PATHS_HEATMAP_OUTPUT_DIR")"
BARPLOT_OUTPUT_DIR="$(normalize_path "$PATHS_BARPLOT_OUTPUT_DIR")"

VISUAL_COLOR_CSV="$(normalize_path "$RESOURCES_VISUAL_COLOR_CSV")"
PCA_VIZ_COLOR_CSV="$(normalize_path "$RESOURCES_PCA_VIZ_COLOR_CSV")"
MAPPING_EXCEL="$(normalize_path "$RESOURCES_MAPPING_EXCEL")"
PROVINCE_MAPPING="$MAPPING_EXCEL"
PHYLOTREE_INDEX_JSON="$(normalize_path "$RESOURCES_PHYLOTREE_INDEX_JSON")"

FORCE_RUN="${PARAMS_FORCE_RUN}"
SAMPLE_META_SHEET="${PARAMS_SAMPLE_META_SHEET}"
ETHNICITY_MAPPING_SHEET="${PARAMS_ETHNICITY_MAPPING_SHEET}"
COUNTRY_MAPPING_SHEET="${PARAMS_COUNTRY_MAPPING_SHEET}"
PROVINCE_MAPPING_SHEET="${PARAMS_PROVINCE_MAPPING_SHEET}"
HEATMAP_GROUP_COL="${PARAMS_HEATMAP_GROUP_COL}"
MIN_SAMPLE_COUNT="${PARAMS_MIN_SAMPLE_COUNT}"
HAPLOGROUP_COLUMN="${PARAMS_HAPLOGROUP_COLUMN}"
PCA_N_COMPONENTS="${PARAMS_PCA_N_COMPONENTS}"
CLASS_BIG_COLUMN="${PARAMS_CLASS_BIG_COLUMN}"
VISUAL_GROUP_COLUMN="${PARAMS_VISUAL_GROUP_COLUMN}"

split_csv_to_array "${PARAMS_NEW_SAMPLE_SOURCE}" NEW_SAMPLE_SOURCE
split_csv_to_array "${PARAMS_COUNTRY_EXTRA_COLS}" COUNTRY_EXTRA_COLS
split_csv_to_array "${PARAMS_ETHNICITY_EXTRA_COLS}" ETHNICITY_EXTRA_COLS
split_csv_to_array "${PARAMS_SELECTED_COLUMNS}" SELECTED_COLUMNS

# 创建必要目录
mkdir -p "$LOG_DIR" "$OUTPUT_DIR" "$DATA_DIR" "$TMP_DIR" "$PCA_OUTPUT_DIR_GLOBAL" "$PCA_OUTPUT_DIR_CHINA" "$PCA_VIS_OUTPUT_DIR" "$HEATMAP_OUTPUT_DIR" "$BARPLOT_OUTPUT_DIR"

# 主日志文件
MAIN_LOG="${LOG_DIR}/run-pipe.sh.log"

# ============================= 2. 日志函数 ============================= #
write_log_file() {
    local level="$1"
    local message="$2"
    local line="[$level] $(date '+%Y-%m-%d %H:%M:%S') - $message"
    echo "$line" >> "$MAIN_LOG"
}

log_info() {
    write_log_file "INFO" "$1"
    ui_info "$1"
}

log_error() {
    write_log_file "ERROR" "$1"
    ui_error "$1"
}

log_success() {
    write_log_file "SUCCESS" "$1"
    ui_ok "$1"
}

# ============================= 3. 参数验证 ============================= #

# 检查base-dir是否存在
if [ ! -d "$BASE_DIR" ]; then
    log_error "基础目录不存在：$BASE_DIR"
    exit 1
fi

# 确保BASE_DIR是绝对路径
BASE_DIR="$(cd "$BASE_DIR" && pwd)"

# PCA核心脚本路径
PCA_CORE_SCRIPT="${PYTHON_DIR}/pca_core.py"

ui_logo
ui_section "mtDNA Haplogroup Pipe" "样本筛选、单倍群整理与群体结构分析"
ui_kv "Config" "$CONFIG_FILE"
ui_kv "BaseDir" "$BASE_DIR"
ui_kv "Input" "$INPUT_EXCEL"
ui_kv "Output" "$OUTPUT_DIR"
ui_kv "Temp" "$TMP_DIR"
ui_kv "Python" "$PYTHON3"
ui_kv "Rscript" "$RSCRIPT"

log_info "========== 开始样本筛选与单倍群分析流程 =========="
log_info "Python环境：$PYTHON3"
log_info "Python环境前缀：$PYTHON_ENV_PREFIX"
log_info "Conda命令：$CONDA_BIN"
log_info "基础目录：$BASE_DIR"
log_info "输入Excel文件：$INPUT_EXCEL"
log_info "数据目录：$DATA_DIR"
log_info "临时文件目录：$TMP_DIR"
if [ "$FORCE_RUN" = "true" ]; then
    log_info "模式：强制重新运行所有步骤"
else
    log_info "模式：跳过已完成步骤"
fi

# ============================= 4. 统一的输入文件检查函数 ============================= #
check_required_files() {
    local missing=0
    
    # 只检查用户需要提供的文件路径
    # conf/ 目录中的配置文件是固定的，在项目部署时已确保存在
    
    # 检查输入Excel文件（用户在base-dir中提供）
    if [ ! -f "$INPUT_EXCEL" ]; then
        log_error "缺失输入文件：$INPUT_EXCEL"
        missing=$((missing + 1))
    fi

    for required_file in \
        "$MAPPING_EXCEL" \
        "$VISUAL_COLOR_CSV" \
        "$PCA_VIZ_COLOR_CSV" \
        "$PHYLOTREE_INDEX_JSON"; do
        if [ ! -f "$required_file" ]; then
            log_error "缺失配置资源文件：$required_file"
            missing=$((missing + 1))
        fi
    done

    if [ ! -x "$PYTHON3" ]; then
        log_error "Python解释器不可执行：$PYTHON3"
        missing=$((missing + 1))
    fi

    if [ ! -x "$RSCRIPT" ]; then
        log_error "Rscript不可执行：$RSCRIPT"
        missing=$((missing + 1))
    fi
    
    if [ $missing -gt 0 ]; then
        log_error "缺失 $missing 个必要文件，无法继续"
        exit 1
    fi
    
    log_info "✅ 所有必要输入文件检查通过"
}

# ============================= 5. 运行状态检查函数 ============================= #
is_step_completed() {
    local step_num=$1
    local step_name=$2
    
    if [ "$FORCE_RUN" = "true" ]; then
        return 1  # 强制运行模式下总是返回未完成
    fi
    
    if grep -q "STEP_${step_num}_COMPLETED" "$MAIN_LOG" 2>/dev/null; then
        return 0  # 已完成
    fi
    return 1  # 未完成
}

mark_step_completed() {
    local step_num=$1
    local step_name=$2
    echo "STEP_${step_num}_COMPLETED" >> "$MAIN_LOG"
    log_success "✅ 完成步骤 $step_num：$step_name"
}

# ============================= 6. 步骤执行函数 ============================= #

# 第一步：单倍群整理
run_step_1() {
    local step_num=1
    local step_name="单倍群整理"
    
    if is_step_completed $step_num "$step_name"; then
        log_info "⊘ 步骤 $step_num 已完成，跳过"
        return 0
    fi
    
    log_info "▶ 开始步骤 $step_num：$step_name"
    
    local PYTHON_SCRIPT="${PYTHON_DIR}/haplogroup_organizer.py"
    
    if ! $PYTHON3 "$PYTHON_SCRIPT" \
        --input "$INPUT_EXCEL" \
        --temp-dir "$TMP_DIR" \
        --data-dir "$OUTPUT_DIR" \
        --tree-index "$PHYLOTREE_INDEX_JSON"; then
        log_error "❌ 步骤 $step_num 执行失败"
        return 1
    fi
    
    mark_step_completed $step_num "$step_name"
}

# 第二步：样本筛选和清洗
run_step_2() {
    local step_num=2
    local step_name="样本筛选和清洗"
    
    if is_step_completed $step_num "$step_name"; then
        log_info "⊘ 步骤 $step_num 已完成，跳过"
        return 0
    fi
    
    log_info "▶ 开始步骤 $step_num：$step_name"
    
    local PYTHON_SCRIPT="${PYTHON_DIR}/data_cleaner.py"
    local OUTPUT_CSV="${OUTPUT_DIR}/本次研究样本基本信息.csv"
    # 检查基础信息速查表
    if [ ! -f "$MAPPING_EXCEL" ]; then
        log_error "缺失基础信息速查表：$MAPPING_EXCEL"
        return 1
    fi
    
    # 分别处理每个NEW_SAMPLE_SOURCE项，合并输出
    local TMP_MERGE="$TMP_DIR/all_new_samples.csv"
    rm -f "$TMP_MERGE"
    local first=1
    for src in "${NEW_SAMPLE_SOURCE[@]}"; do
        local TMP_OUT="$TMP_DIR/new_sample_${src}.csv"
        if ! $PYTHON3 "$PYTHON_SCRIPT" \
            --input "$INPUT_EXCEL" \
            --mapping "$MAPPING_EXCEL" \
            --new-source "$src" \
            --sample-meta-sheet "$SAMPLE_META_SHEET" \
            --ethnicity-mapping-sheet "$ETHNICITY_MAPPING_SHEET" \
            --country-mapping-sheet "$COUNTRY_MAPPING_SHEET" \
            --province-mapping-sheet "$PROVINCE_MAPPING_SHEET" \
            --country-cols "${COUNTRY_EXTRA_COLS[@]}" \
            --ethnicity-cols "${ETHNICITY_EXTRA_COLS[@]}" \
            --output "$TMP_OUT"; then
            log_error "❌ 步骤 $step_num 执行失败（source=$src）"
            return 1
        fi
        if [ $first -eq 1 ]; then
            cat "$TMP_OUT" > "$TMP_MERGE"
            first=0
        else
            tail -n +2 "$TMP_OUT" >> "$TMP_MERGE"
        fi
    done
    mkdir -p "$OUTPUT_DIR"
    if ! mv "$TMP_MERGE" "$OUTPUT_CSV"; then
        log_error "❌ 步骤 $step_num 执行失败：无法写入输出文件 $OUTPUT_CSV"
        return 1
    fi
    mark_step_completed $step_num "$step_name"
}

# 第三步：合并单倍群信息
run_step_3() {
    local step_num=3
    local step_name="合并单倍群信息"
    
    if is_step_completed $step_num "$step_name"; then
        log_info "⊘ 步骤 $step_num 已完成，跳过"
        return 0
    fi
    
    log_info "▶ 开始步骤 $step_num：$step_name"
    
    local PYTHON_SCRIPT="${PYTHON_DIR}/merger.py"
    local HAPLOGROUP_CSV="${OUTPUT_DIR}/全部单倍群整理.csv"
    local SAMPLE_CSV="${OUTPUT_DIR}/本次研究样本基本信息.csv"
    local OUTPUT_CSV="${OUTPUT_DIR}/本次研究样本基本信息_含单倍群.csv"
    
    # 检查输入文件
    if [ ! -f "$HAPLOGROUP_CSV" ]; then
        log_error "缺失单倍群文件（需先完成第1步）：$HAPLOGROUP_CSV"
        return 1
    fi
    
    if [ ! -f "$SAMPLE_CSV" ]; then
        log_error "缺失样本信息文件（需先完成第2步）：$SAMPLE_CSV"
        return 1
    fi
    
    if ! $PYTHON3 "$PYTHON_SCRIPT" \
        --haplogroup "$HAPLOGROUP_CSV" \
        --sample "$SAMPLE_CSV" \
        --output "$OUTPUT_CSV"; then
        log_error "❌ 步骤 $step_num 执行失败"
        return 1
    fi
    
    mark_step_completed $step_num "$step_name"
}

# 第四步：频率分析前置处理（含中国数据分离）
run_step_4() {
    local step_num=4
    local step_name="频率分析前置处理（含中国数据分离）"
    
    if is_step_completed $step_num "$step_name"; then
        log_info "⊘ 步骤 $step_num 已完成，跳过"
        return 0
    fi
    
    log_info "▶ 开始步骤 $step_num：$step_name"
    
    local PYTHON_SCRIPT="${PYTHON_DIR}/frequency_preprocessor.py"
    local INPUT_CSV="${OUTPUT_DIR}/本次研究样本基本信息_含单倍群.csv"
    local OUTPUT_PREPARED="${OUTPUT_DIR}/频率分析准备.csv"
    local OUTPUT_PREPARED_CHINA="${OUTPUT_DIR}/频率分析准备_中国.csv"
    local OUTPUT_FREQUENCY="${OUTPUT_DIR}/频率矩阵.csv"
    local OUTPUT_FREQUENCY_CHINA="${OUTPUT_DIR}/频率矩阵_中国.csv"
    
    # 检查输入文件
    if [ ! -f "$INPUT_CSV" ]; then
        log_error "缺失输入文件（需先完成第3步）：$INPUT_CSV"
        return 1
    fi
    
    if [ ! -f "$PROVINCE_MAPPING" ]; then
        log_error "缺失省份映射表：$PROVINCE_MAPPING"
        return 1
    fi
    
    if ! $PYTHON3 "$PYTHON_SCRIPT" \
        --input "$INPUT_CSV" \
        --province-mapping "$PROVINCE_MAPPING" \
        --selected-cols "${SELECTED_COLUMNS[@]}" \
        --min-count "$MIN_SAMPLE_COUNT" \
        --haplogroup-col "$HAPLOGROUP_COLUMN" \
        --output-prepared "$OUTPUT_PREPARED" \
        --output-prepared-china "$OUTPUT_PREPARED_CHINA" \
        --output-frequency "$OUTPUT_FREQUENCY" \
        --output-frequency-china "$OUTPUT_FREQUENCY_CHINA"; then
        log_error "❌ 步骤 $step_num 执行失败"
        return 1
    fi

    # 补充：在生成频率准备文件后，直接提取热图分组映射（分别对全球与中国文件）
    local HM_SCRIPT="${PYTHON_DIR}/heatmap_group_extractor.py"
    local HM_INPUT_PREPARED="$OUTPUT_PREPARED"

    if [ ! -f "$HM_INPUT_PREPARED" ]; then
        log_error "缺失准备数据（已完成第4步但未找到输出）：$HM_INPUT_PREPARED"
        return 1
    fi

    mkdir -p "$HEATMAP_OUTPUT_DIR"

        # 提取全球映射
        if ! $PYTHON3 "$HM_SCRIPT" \
            --input "$HM_INPUT_PREPARED" \
            --group-col "$HEATMAP_GROUP_COL" \
            --output-dir "$HEATMAP_OUTPUT_DIR"; then
        log_error "❌ 步骤 $step_num 补充任务失败（热图分组提取）"
        return 1
    fi

    # 提取中国映射（基于中国准备文件，不使用筛选）
    local HM_INPUT_PREPARED_CHINA="$OUTPUT_PREPARED_CHINA"
    if [ -f "$HM_INPUT_PREPARED_CHINA" ]; then
        if ! $PYTHON3 "$HM_SCRIPT" \
            --input "$HM_INPUT_PREPARED_CHINA" \
            --group-col "$HEATMAP_GROUP_COL" \
            --suffix "_China" \
            --output-dir "$HEATMAP_OUTPUT_DIR"; then
            log_error "❌ 步骤 $step_num 补充任务失败（热图分组提取，中国）"
            return 1
        fi
    else
        log_error "缺失中国准备数据（需启用中国输出）：$HM_INPUT_PREPARED_CHINA"
        return 1
    fi
    
    mark_step_completed $step_num "$step_name"
}

# 第五步：PCA降维（全球）
run_step_5() {
    local step_num=5
    local step_name="PCA降维（全球）"
    
    if is_step_completed $step_num "$step_name"; then
        log_info "⊘ 步骤 $step_num 已完成，跳过"
        return 0
    fi
    
    log_info "▶ 开始步骤 $step_num：$step_name"
    
    local FREQUENCY_CSV="${OUTPUT_DIR}/频率矩阵.csv"
    
    # 检查输入文件
    if [ ! -f "$FREQUENCY_CSV" ]; then
        log_error "缺失频率表（需先完成第4步）：$FREQUENCY_CSV"
        return 1
    fi
    
    if [ ! -f "$PCA_CORE_SCRIPT" ]; then
        log_error "缺失PCA核心脚本：$PCA_CORE_SCRIPT"
        return 1
    fi
    
    # 创建PCA输出目录
    mkdir -p "$PCA_OUTPUT_DIR_GLOBAL"
    
    # 直接调用PCA核心脚本
    if ! $PYTHON3 "$PCA_CORE_SCRIPT" \
        --input "$FREQUENCY_CSV" \
        --output-dir "$PCA_OUTPUT_DIR_GLOBAL" \
        --n-components "$PCA_N_COMPONENTS"; then
        log_error "❌ 步骤 $step_num 执行失败"
        return 1
    fi
    
    # 验证PCA结果是否生成（不再复制到 data 目录）
    if [ -f "$PCA_OUTPUT_DIR_GLOBAL/pca_results_${PCA_N_COMPONENTS}.csv" ]; then
        :
    else
        log_error "❌ PCA输出文件不存在"
        return 1
    fi
    
    mark_step_completed $step_num "$step_name"
}

# 第六步：PCA后处理（全球）
run_step_6() {
    local step_num=6
    local step_name="PCA后处理（全球）"
    
    if is_step_completed $step_num "$step_name"; then
        log_info "⊘ 步骤 $step_num 已完成，跳过"
        return 0
    fi
    
    log_info "▶ 开始步骤 $step_num：$step_name"
    
    local PYTHON_SCRIPT="${PYTHON_DIR}/pca_postprocessor.py"
    local PCA_RESULT_FILE="${PCA_OUTPUT_DIR_GLOBAL}/pca_results_${PCA_N_COMPONENTS}.csv"
    local PREPARED_CSV="${OUTPUT_DIR}/频率分析准备.csv"
    local PCA_VIS_CSV="${PCA_OUTPUT_DIR_GLOBAL}/PCA输入文件.csv"
    
    # 检查输入文件
    if [ ! -f "$PCA_RESULT_FILE" ]; then
        log_error "缺失PCA结果（需先完成第5步）：$PCA_RESULT_FILE"
        return 1
    fi
    
    if [ ! -f "$PREPARED_CSV" ]; then
        log_error "缺失准备数据（需先完成第4步）：$PREPARED_CSV"
        return 1
    fi
    
    if ! $PYTHON3 "$PYTHON_SCRIPT" \
        --pca-result "$PCA_RESULT_FILE" \
        --prepared-data "$PREPARED_CSV" \
        --metadata-cols "Classification" "Continent" "Group" \
        --group-column "Group" \
        --class-big-column "$CLASS_BIG_COLUMN" \
        --output "$PCA_VIS_CSV"; then
        log_error "❌ 步骤 $step_num 执行失败"
        return 1
    fi
    
    mark_step_completed $step_num "$step_name"
}

# 第七步：PCA降维（中国）
run_step_7() {
    local step_num=7
    local step_name="PCA降维（中国）"
    
    if is_step_completed $step_num "$step_name"; then
        log_info "⊘ 步骤 $step_num 已完成，跳过"
        return 0
    fi
    
    log_info "▶ 开始步骤 $step_num：$step_name"
    
    local CHINA_FREQUENCY="${OUTPUT_DIR}/频率矩阵_中国.csv"
    
    # 检查输入文件
    if [ ! -f "$CHINA_FREQUENCY" ]; then
        log_error "缺失中国频率表（需先完成第4步）：$CHINA_FREQUENCY"
        return 1
    fi
    
    if [ ! -f "$PCA_CORE_SCRIPT" ]; then
        log_error "缺失PCA核心脚本：$PCA_CORE_SCRIPT"
        return 1
    fi
    
    # 创建PCA输出目录
    mkdir -p "$PCA_OUTPUT_DIR_CHINA"
    
    # 直接调用PCA核心脚本
    if ! $PYTHON3 "$PCA_CORE_SCRIPT" \
        --input "$CHINA_FREQUENCY" \
        --output-dir "$PCA_OUTPUT_DIR_CHINA" \
        --n-components "$PCA_N_COMPONENTS"; then
        log_error "❌ 步骤 $step_num 执行失败"
        return 1
    fi
    
    # 验证PCA结果是否生成（不再复制到 data 目录）
    if [ -f "$PCA_OUTPUT_DIR_CHINA/pca_results_${PCA_N_COMPONENTS}.csv" ]; then
        :
    else
        log_error "❌ PCA输出文件不存在"
        return 1
    fi
    
    mark_step_completed $step_num "$step_name"
}

# 第八步：PCA后处理（中国）
run_step_8() {
    local step_num=8
    local step_name="PCA后处理（中国）"
    
    if is_step_completed $step_num "$step_name"; then
        log_info "⊘ 步骤 $step_num 已完成，跳过"
        return 0
    fi
    
    log_info "▶ 开始步骤 $step_num：$step_name"
    
    local PYTHON_SCRIPT="${PYTHON_DIR}/pca_postprocessor.py"
    local CHINA_PCA_RESULT_FILE="${PCA_OUTPUT_DIR_CHINA}/pca_results_${PCA_N_COMPONENTS}.csv"
    local CHINA_PREPARED="${OUTPUT_DIR}/频率分析准备_中国.csv"
    local CHINA_PCA_VIS="${PCA_OUTPUT_DIR_CHINA}/PCA输入文件.csv"
    
    # 检查输入文件
    if [ ! -f "$CHINA_PCA_RESULT_FILE" ]; then
        log_error "缺失中国PCA结果（需先完成第7步）：$CHINA_PCA_RESULT_FILE"
        return 1
    fi
    
    if [ ! -f "$CHINA_PREPARED" ]; then
        log_error "缺失中国准备数据（需先完成第4步）：$CHINA_PREPARED"
        return 1
    fi
    
    if ! $PYTHON3 "$PYTHON_SCRIPT" \
        --pca-result "$CHINA_PCA_RESULT_FILE" \
        --prepared-data "$CHINA_PREPARED" \
        --metadata-cols "Classification" "Group" \
        --group-column "Group" \
        --class-big-column "$CLASS_BIG_COLUMN" \
        --output "$CHINA_PCA_VIS"; then
        log_error "❌ 步骤 $step_num 执行失败"
        return 1
    fi
    
    mark_step_completed $step_num "$step_name"
}

# 第九步：频率柱状图（全球）
run_step_9() {
    local step_num=9
    local step_name="频率柱状图（全球+中国）"

    if is_step_completed $step_num "$step_name"; then
        log_info "⊘ 步骤 $step_num 已完成，跳过"
        return 0
    fi

    log_info "▶ 开始步骤 $step_num：$step_name"

    local PYTHON_SCRIPT="${PYTHON_DIR}/haplogroup_barplot.py"

    # 全局输入/输出
    local INPUT_CSV_GLOBAL="${OUTPUT_DIR}/频率矩阵.csv"
    local PREPARED_CSV_GLOBAL="${OUTPUT_DIR}/频率分析准备.csv"
    local OUTPUT_PDF_GLOBAL="${BARPLOT_OUTPUT_DIR}/Haplogroup_distribution.pdf"

    # 中国输入/输出
    local INPUT_CSV_CHINA="${OUTPUT_DIR}/频率矩阵_中国.csv"
    local PREPARED_CSV_CHINA="${OUTPUT_DIR}/频率分析准备_中国.csv"
    local OUTPUT_PDF_CHINA="${BARPLOT_OUTPUT_DIR}/Haplogroup_distribution_中国.pdf"

    # 检查颜色映射文件
    if [ ! -f "$VISUAL_COLOR_CSV" ]; then
        log_error "缺失颜色映射文件：$VISUAL_COLOR_CSV"
        return 1
    fi

    mkdir -p "$BARPLOT_OUTPUT_DIR"

    # 先全球
    if [ ! -f "$INPUT_CSV_GLOBAL" ]; then
        log_error "缺失全球频率表（需先完成第4步）：$INPUT_CSV_GLOBAL"
        return 1
    fi
    if [ ! -f "$PREPARED_CSV_GLOBAL" ]; then
        log_error "缺失全球准备数据（需先完成第4步）：$PREPARED_CSV_GLOBAL"
        return 1
    fi
    if ! $PYTHON3 "$PYTHON_SCRIPT" \
        --input "$INPUT_CSV_GLOBAL" \
        --prepared "$PREPARED_CSV_GLOBAL" \
        --color "$VISUAL_COLOR_CSV" \
        --group-col "$VISUAL_GROUP_COLUMN" \
        --output-pdf "$OUTPUT_PDF_GLOBAL"; then
        log_error "❌ 全球柱状图生成失败"
        return 1
    fi

    # 再中国
    if [ ! -f "$INPUT_CSV_CHINA" ]; then
        log_error "缺失中国频率表（需先完成第4步）：$INPUT_CSV_CHINA"
        return 1
    fi
    if [ ! -f "$PREPARED_CSV_CHINA" ]; then
        log_error "缺失中国准备数据（需先完成第4步并启用中国输出）：$PREPARED_CSV_CHINA"
        return 1
    fi
    if ! $PYTHON3 "$PYTHON_SCRIPT" \
        --input "$INPUT_CSV_CHINA" \
        --prepared "$PREPARED_CSV_CHINA" \
        --color "$VISUAL_COLOR_CSV" \
        --group-col "$VISUAL_GROUP_COLUMN" \
        --output-pdf "$OUTPUT_PDF_CHINA"; then
        log_error "❌ 中国柱状图生成失败"
        return 1
    fi

    mark_step_completed $step_num "$step_name"
}

#（已移除）原第十步：提取热图分组映射

# 第十步：PCA可视化（R脚本，全球+中国）
run_step_10() {
    local step_num=10
    local step_name="PCA可视化（R脚本，全球+中国）"

    if is_step_completed $step_num "$step_name"; then
        log_info "⊘ 步骤 $step_num 已完成，跳过"
        return 0
    fi

    log_info "▶ 开始步骤 $step_num：$step_name"

    local R_SCRIPT="${SRC_DIR}/pca_visualization.R"
    local PCA_OUTPUT_DIR="${PCA_VIS_OUTPUT_DIR}"

    if ! command -v "$RSCRIPT" &> /dev/null; then
        log_error "缺失Rscript命令，请确认R环境已安装"
        return 1
    fi

    if [ ! -f "$R_SCRIPT" ]; then
        log_error "缺失R脚本：$R_SCRIPT"
        return 1
    fi

    mkdir -p "$PCA_OUTPUT_DIR"

    # 全球PCA可视化
    local PCA_INPUT_GLOBAL="${PCA_OUTPUT_DIR_GLOBAL}/PCA输入文件.csv"
    
    if [ ! -f "$PCA_INPUT_GLOBAL" ]; then
        log_error "缺失全球PCA输入（需先完成第6步）：$PCA_INPUT_GLOBAL"
        return 1
    fi
    
    if [ ! -f "$PCA_VIZ_COLOR_CSV" ]; then
        log_error "缺失PCA颜色映射文件：$PCA_VIZ_COLOR_CSV"
        return 1
    fi

    log_info "  运行全球PCA可视化..."
    if ! "$RSCRIPT" "$R_SCRIPT" \
        --pca-input "$PCA_INPUT_GLOBAL" \
        --color-file "$PCA_VIZ_COLOR_CSV" \
        --output-dir "$PCA_OUTPUT_DIR" \
        --output-name "PCA_visualization_global.pdf"; then
        log_error "❌ 全球PCA可视化失败"
        return 1
    fi

    # 中国PCA可视化
    local PCA_INPUT_CHINA="${PCA_OUTPUT_DIR_CHINA}/PCA输入文件.csv"
    
    if [ ! -f "$PCA_INPUT_CHINA" ]; then
        log_error "缺失中国PCA输入（需先完成第8步）：$PCA_INPUT_CHINA"
        return 1
    fi

    log_info "  运行中国PCA可视化..."
    if ! "$RSCRIPT" "$R_SCRIPT" \
        --pca-input "$PCA_INPUT_CHINA" \
        --color-file "$PCA_VIZ_COLOR_CSV" \
        --output-dir "$PCA_OUTPUT_DIR" \
        --output-name "PCA_visualization_china.pdf"; then
        log_error "❌ 中国PCA可视化失败"
        return 1
    fi

    mark_step_completed $step_num "$step_name"
}

# 第十一步：频率矩阵热图可视化（R脚本，全球+中国）
run_step_11() {
    local step_num=11
    local step_name="频率矩阵热图可视化（R脚本，全球+中国）"

    if is_step_completed $step_num "$step_name"; then
        log_info "⊘ 步骤 $step_num 已完成，跳过"
        return 0
    fi

    log_info "▶ 开始步骤 $step_num：$step_name"

    local R_SCRIPT="${SRC_DIR}/heatmap_visualization.R"
    local HEATMAP_OUTPUT_DIR="${PATHS_HEATMAP_OUTPUT_DIR}"
    HEATMAP_OUTPUT_DIR="$(normalize_path "$HEATMAP_OUTPUT_DIR")"

    if ! command -v "$RSCRIPT" &> /dev/null; then
        log_error "缺失Rscript命令，请确认R环境已安装"
        return 1
    fi

    if [ ! -f "$R_SCRIPT" ]; then
        log_error "缺失R脚本：$R_SCRIPT"
        return 1
    fi

    mkdir -p "$HEATMAP_OUTPUT_DIR"

    local FREQ_MATRIX_GLOBAL="${OUTPUT_DIR}/频率矩阵.csv"
    local GROUP_MAPPING_GLOBAL="${HEATMAP_OUTPUT_DIR}/Classification_Group_mapping.csv"
    
    if [ ! -f "$FREQ_MATRIX_GLOBAL" ]; then
        log_error "缺失全球频率矩阵（需先完成第4步）：$FREQ_MATRIX_GLOBAL"
        return 1
    fi
    
    if [ ! -f "$GROUP_MAPPING_GLOBAL" ]; then
        log_error "缺失全球分组映射（需先完成第4步）：$GROUP_MAPPING_GLOBAL"
        return 1
    fi

    log_info "  运行全球频率热图..."
    if ! "$RSCRIPT" "$R_SCRIPT" \
        --frequency-matrix "$FREQ_MATRIX_GLOBAL" \
        --group-mapping "$GROUP_MAPPING_GLOBAL" \
        --output-dir "$HEATMAP_OUTPUT_DIR" \
        --output-prefix "heatmap_global"; then
        log_error "❌ 全球热图可视化失败"
        return 1
    fi

    local FREQ_MATRIX_CHINA="${OUTPUT_DIR}/频率矩阵_中国.csv"
    local GROUP_MAPPING_CHINA="${HEATMAP_OUTPUT_DIR}/Classification_Group_mapping_China.csv"
    
    if [ ! -f "$FREQ_MATRIX_CHINA" ]; then
        log_error "缺失中国频率矩阵（需先完成第4步）：$FREQ_MATRIX_CHINA"
        return 1
    fi
    
    if [ ! -f "$GROUP_MAPPING_CHINA" ]; then
        log_error "缺失中国分组映射（需先完成第4步）：$GROUP_MAPPING_CHINA"
        return 1
    fi

    log_info "  运行中国频率热图..."
    if ! "$RSCRIPT" "$R_SCRIPT" \
        --frequency-matrix "$FREQ_MATRIX_CHINA" \
        --group-mapping "$GROUP_MAPPING_CHINA" \
        --output-dir "$HEATMAP_OUTPUT_DIR" \
        --output-prefix "heatmap_china"; then
        log_error "❌ 中国热图可视化失败"
        return 1
    fi

    mark_step_completed $step_num "$step_name"
}

main() {
    # 检查所有输入文件
    check_required_files
    
    # 执行步骤
    run_step_1 || exit 1
    run_step_2 || exit 1
    run_step_3 || exit 1
    run_step_4 || exit 1  # 频率处理（同时输出全球和中国）
    run_step_5 || exit 1  # 全球PCA
    run_step_6 || exit 1  # 全球PCA后处理
    run_step_7 || exit 1  # 中国PCA
    run_step_8 || exit 1  # 中国PCA后处理
    run_step_9 || exit 1   # 频率柱状图（全球+中国）
    run_step_10 || exit 1  # PCA可视化（全球+中国）
    run_step_11 || exit 1  # 频率热图可视化（全球+中国）
    
    log_success "========== 流程执行完成 =========="
}

# 运行主流程
main
