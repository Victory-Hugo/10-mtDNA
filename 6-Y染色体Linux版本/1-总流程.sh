#!/bin/bash

################################################################################
# 脚本名: run_haplo_app.sh
# 功能: 单倍群频率分析 - Haplogroup Population Profile Generator
# 作者: BigLin
# 创建日期: 2025-10-31
# 描述: 将R脚本转换为Linux Shell脚本，支持并行处理和断点续跑
################################################################################

set -o pipefail

################################################################################
# 第一部分: 依赖检查与环境设置
################################################################################

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 【用户配置区】所有可配置的变量都在这里，请按需修改
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# 输入/输出路径配置 - 全部写死，请修改为实际路径
readonly INPUT_FILE="/mnt/c/Users/Administrator/Desktop/Run_HaplogroupPPG_APP/input/Table_Haplogroup.xlsx"
readonly EXCEL_SHEET="MergeOneHanwithAncients"
readonly OUTPUT_DIR="/mnt/c/Users/Administrator/Desktop/Run_HaplogroupPPG_APP/output"

# NC 范围配置 
readonly NC_ARRAY=(1 2 3 4 5 6 7 8)

# 处理参数配置
readonly USE_EQUALIZE=false
readonly THR_COMMON=0.05
readonly THR_RARE=0.01
readonly JSD_EPS=1e-12

# 并行处理配置
readonly MAX_PARALLEL=4  # 最多同时运行4个NC分析

# 绘图参数配置
readonly PLOT_BASE_WIDTH=8
readonly PLOT_BASE_HEIGHT=6

# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
# 【系统要求】一般不需要修改
# ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━

# 前置要求
REQUIRED_TOOLS=("Rscript" "awk" "cut" "paste" "sort" "uniq" "parallel" "tput")
REQUIRED_R_PACKAGES=("readxl" "dplyr" "tidyr" "stringr" "ggplot2" "ggrepel" "FactoMineR" "vegan" "pheatmap" "ggalluvial" "scales" "ape" "grid")

cat << 'EOF'
╔════════════════════════════════════════════════════════════════════════════╗
║              Haplogroup Population Profile Generator v2.0                  ║
║                           Linux Shell 版本                                 ║
║                         Author: BigLin (2025)                              ║
╚════════════════════════════════════════════════════════════════════════════╝

EOF

# 颜色定义
readonly RED='\033[0;31m'
readonly GREEN='\033[0;32m'
readonly YELLOW='\033[1;33m'
readonly BLUE='\033[0;34m'
readonly CYAN='\033[0;36m'
readonly BOLD='\033[1m'
readonly NC='\033[0m' # No Color

# 日志函数
log_info() {
    echo -e "${BLUE}[$(date '+%H:%M:%S')]${NC} ${GREEN}✓${NC} $1"
}

log_warn() {
    echo -e "${YELLOW}[$(date '+%H:%M:%S')] ⚠ WARNING:${NC} $1"
}

log_error() {
    echo -e "${RED}[$(date '+%H:%M:%S')] ✗ ERROR:${NC} $1"
}

log_section() {
    echo -e "${BOLD}${CYAN}════════════════════════════════════════════════════════════${NC}"
    echo -e "${BOLD}${CYAN}► $1${NC}"
    echo -e "${BOLD}${CYAN}════════════════════════════════════════════════════════════${NC}"
}

# 进度条函数
show_progress() {
    local current=$1
    local total=$2
    local task=$3
    local percent=$((current * 100 / total))
    local filled=$((percent / 5))
    local empty=$((20 - filled))
    
    printf "\r${BLUE}[%-20s]${NC} ${YELLOW}%3d%%${NC} (%d/%d) ${GREEN}%s${NC}" \
        "$(printf '#%.0s' $(seq 1 $filled))$(printf '·%.0s' $(seq 1 $empty))" \
        "$percent" "$current" "$total" "$task"
}

# 信号处理 - 清理未完成的任务
cleanup() {
    log_error "脚本被中断！"
    
    if [ -n "$TEMP_DIR" ] && [ -d "$TEMP_DIR" ]; then
        log_warn "清理临时文件..."
        rm -rf "$TEMP_DIR"
    fi
    
    # 终止所有后台进程
    if [ -n "$PID_LOG" ] && [ -f "$PID_LOG" ]; then
        while read -r pid; do
            if kill -0 "$pid" 2>/dev/null; then
                kill -9 "$pid" 2>/dev/null
            fi
        done < "$PID_LOG"
        rm -f "$PID_LOG"
    fi
    
    log_warn "已保留完成的任务，已清理未完成的任务"
    exit 130
}

trap cleanup SIGINT SIGTERM

################################################################################
# 第二部分: 依赖检查
################################################################################

check_dependencies() {
    log_section "检查系统依赖"
    
    local missing_tools=()
    
    for tool in "${REQUIRED_TOOLS[@]}"; do
        if ! command -v "$tool" &> /dev/null; then
            missing_tools+=("$tool")
            log_warn "缺少: $tool"
        else
            log_info "已找到: $tool"
        fi
    done
    
    if [ ${#missing_tools[@]} -gt 0 ]; then
        log_error "缺少以下工具: ${missing_tools[*]}"
        log_error "请在 Linux 上运行: sudo apt-get install parallel moreutils"
        exit 1
    fi
    
    # 检查R和R包
    log_section "检查 R 环境"
    
    if ! command -v Rscript &> /dev/null; then
        log_error "找不到 Rscript！请先安装 R"
        exit 1
    fi
    
    R_VERSION=$(Rscript --version 2>&1 | head -1)
    log_info "检测到: $R_VERSION"
    
    # 检查R包
    log_info "检查 R 包依赖..."
    
    R_INSTALL_SCRIPT="
library_check <- function() {
    pkgs <- c('readxl','dplyr','tidyr','stringr','ggplot2','ggrepel','FactoMineR','vegan','pheatmap','ggalluvial','scales','ape','grid')
    missing <- c()
    for (p in pkgs) {
        if (!requireNamespace(p, quietly=TRUE)) {
            missing <- c(missing, p)
        }
    }
    if (length(missing) > 0) {
        cat('MISSING:', paste(missing, collapse=' '), '\n')
    } else {
        cat('OK\n')
    }
}
library_check()
"
    
    MISSING_PKGS=$(Rscript -e "$R_INSTALL_SCRIPT" 2>/dev/null | grep "MISSING" | sed 's/MISSING: //')
    
    if [ -n "$MISSING_PKGS" ]; then
        log_warn "缺少 R 包: $MISSING_PKGS"
        log_info "自动安装缺失的 R 包..."
        
        for pkg in $MISSING_PKGS; do
            Rscript -e "install.packages('$pkg', repos='https://cloud.r-project.org')" 2>&1 | tail -1
        done
    else
        log_info "所有 R 包已安装"
    fi
}

################################################################################
# 第三部分: 环境初始化（基于第一部分的配置）
################################################################################

configure_variables() {
    log_section "初始化分析环境"
    
    # 获取脚本所在目录并切换到该目录
    readonly SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
    readonly WORK_DIR="$SCRIPT_DIR"
    
    # 切换到脚本目录（确保R脚本可以被找到）
    cd "$WORK_DIR" || {
        log_error "无法切换到脚本目录: $WORK_DIR"
        exit 1
    }
    
    # 并行处理配置
    readonly TEMP_DIR="/tmp/haplo_$$"
    readonly PID_LOG="$TEMP_DIR/pids.log"
    readonly PROGRESS_LOG="$TEMP_DIR/progress.log"
    readonly STATE_DIR="$TEMP_DIR/states"
    
    mkdir -p "$TEMP_DIR" "$STATE_DIR"
    touch "$PID_LOG" "$PROGRESS_LOG"
    
    log_info "脚本目录: $SCRIPT_DIR"
    log_info "工作目录: $WORK_DIR"
    log_info "输入文件: $INPUT_FILE"
    log_info "Excel工作表: $EXCEL_SHEET"
    log_info "输出目录: $OUTPUT_DIR"
    log_info "并行数: $MAX_PARALLEL"
    log_info "NC 范围: ${NC_ARRAY[*]}"
}

################################################################################
# 第四部分: 输入文件验证
################################################################################

verify_input_file() {
    local file_path="$1"
    
    if [ ! -f "$file_path" ]; then
        log_error "输入文件不存在: $file_path"
        return 1
    fi
    
    log_info "输入文件已验证: $file_path"
    return 0
}

################################################################################
# 第五部分: 核心分析流程管理
################################################################################

# 运行单个NC的完整分析流程
run_nc_analysis() {
    local input_file="$1"
    local nc="$2"
    local output_dir="$3"
    local temp_dir="$4"
    
    # 创建NC专用的输出目录
    local nc_output_dir="$output_dir/out_NC$(printf '%02d' $nc)"
    mkdir -p "$nc_output_dir"
    
    # 步骤1: 数据预处理
    local processed_data="$temp_dir/processed_data_NC${nc}.rds"
    if ! Rscript "2-data_preprocessing.R" "$input_file" "$processed_data" "$nc" "$USE_EQUALIZE" "$EXCEL_SHEET" 2>&1; then
        log_error "数据预处理失败 (NC=$nc)"
        return 1
    fi
    
    # 步骤2: 频率计算
    if ! Rscript "3-frequency_calculation.R" "$processed_data" "$nc_output_dir" "$nc" 2>&1; then
        log_error "频率计算失败 (NC=$nc)"
        return 1
    fi
    
    # 步骤3: 距离计算
    local freq_matrix="$nc_output_dir/tables/freq_matrix_NC$(printf '%02d' $nc).csv"
    if ! Rscript "4-distance_calculation.R" "$freq_matrix" "$nc_output_dir" "$nc" 2>&1; then
        log_error "距离计算失败 (NC=$nc)"
        return 1
    fi
    
    # 步骤4: 可视化
    local distance_matrices="$nc_output_dir/tables/distance_matrices_NC$(printf '%02d' $nc).rds"
    if ! Rscript "5-visualization.R" "$distance_matrices" "$freq_matrix" "$nc_output_dir" "$nc" 2>&1; then
        log_error "可视化失败 (NC=$nc)"
        return 1
    fi
    
    return 0
}

################################################################################
# 第六部分: 并行处理管理
################################################################################

process_nc_parallel() {
    local nc_value=$1
    local input_file=$2
    local state_file="$STATE_DIR/NC${nc_value}.state"
    
    # 检查是否已完成
    if [ -f "$state_file" ]; then
        local state=$(cat "$state_file")
        if [ "$state" = "COMPLETED" ]; then
            log_info "NC=$nc_value 已完成（断点续跑）"
            return 0
        fi
    fi
    
    # 运行分析流程
    {
        if run_nc_analysis "$input_file" "$nc_value" "$OUTPUT_DIR" "$TEMP_DIR" 2>&1; then
            echo "COMPLETED" > "$state_file"
            log_info "✓ NC=$nc_value 处理完成"
        else
            echo "FAILED" > "$state_file"
            log_error "✗ NC=$nc_value 处理失败"
            return 1
        fi
    }
}

run_parallel_analysis() {
    local input_file="$1"
    log_section "并行处理 NC 分析"
    
    local total=${#NC_ARRAY[@]}
    local completed=0
    local failed=0
    local pids=()
    
    # 启动所有后台进程
    for nc in "${NC_ARRAY[@]}"; do
        process_nc_parallel "$nc" "$input_file" &
        local pid=$!
        pids+=($pid)
        echo $pid >> "$PID_LOG"
        
        show_progress $((completed + 1)) $total "启动 NC=$nc (PID=$pid)"
        ((completed++))
        
        # 限制并发数
        if [ $((completed % MAX_PARALLEL)) -eq 0 ]; then
            sleep 1
        fi
    done
    
    # 等待所有进程完成
    completed=0
    for pid in "${pids[@]}"; do
        if wait "$pid"; then
            ((completed++))
        else
            ((failed++))
        fi
        show_progress $((completed + failed)) $total "等待完成..."
    done
    
    echo "" # 换行
    
    if [ $failed -eq 0 ]; then
        log_info "所有 NC 分析已完成！"
        return 0
    else
        log_warn "有 $failed 个任务失败"
        return 1
    fi
}

################################################################################
# 第七部分: 后处理与汇总
################################################################################

post_processing() {
    log_section "后处理与汇总"
    
    if [ ! -d "$OUTPUT_DIR" ]; then
        log_error "输出目录不存在"
        return 1
    fi
    
    # 生成汇总文件
    local summary_file="$OUTPUT_DIR/ANALYSIS_SUMMARY.txt"
    
    {
        echo "╔════════════════════════════════════════════════════════════╗"
        echo "║         单倍群频率分析 - 汇总报告                           ║"
        echo "╚════════════════════════════════════════════════════════════╝"
        echo ""
        echo "分析时间: $(date '+%Y-%m-%d %H:%M:%S')"
        echo "输入文件: $INPUT_FILE"
        echo "输出目录: $OUTPUT_DIR"
        echo ""
        echo "处理的 NC 范围: ${NC_ARRAY[*]}"
        echo ""
        echo "生成的文件结构:"
        find "$OUTPUT_DIR" -type d | head -20
        echo ""
        echo "详细内容:"
        find "$OUTPUT_DIR" -type f -name "*.csv" | wc -l | xargs echo "CSV 文件数:"
        find "$OUTPUT_DIR" -type f -name "*.pdf" | wc -l | xargs echo "PDF 文件数:"
        find "$OUTPUT_DIR" -type f -name "*.png" | wc -l | xargs echo "PNG 文件数:"
    } > "$summary_file"
    
    log_info "汇总报告已保存至: $summary_file"
    
    # 打印汇总
    cat "$summary_file"
}

################################################################################
# 第八部分: 主程序流程
################################################################################

main() {
    local start_time=$(date +%s)
    
    # 初始化
    configure_variables
    check_dependencies
    
    # 验证输入文件
    log_section "验证输入文件"
    if ! verify_input_file "$INPUT_FILE"; then
        exit 1
    fi
    
    # 并行处理
    run_parallel_analysis "$INPUT_FILE"
    local analysis_result=$?
    
    # 后处理
    post_processing
    
    # 清理临时文件
    log_section "清理临时文件"
    rm -rf "$TEMP_DIR"
    log_info "临时文件已清理"
    
    # 计算运行时间
    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    local hours=$((duration / 3600))
    local minutes=$(((duration % 3600) / 60))
    local seconds=$((duration % 60))
    
    # 最终输出
    log_section "分析完成！"
    
    echo -e "${BOLD}${GREEN}════════════════════════════════════════════════════════════${NC}"
    echo -e "${BOLD}${GREEN}✓ 所有分析已完成！${NC}"
    echo -e "${BOLD}${GREEN}════════════════════════════════════════════════════════════${NC}"
    echo ""
    echo "运行时长: ${hours}h ${minutes}m ${seconds}s"
    echo "结果位置: $(cd "$OUTPUT_DIR" && pwd)"
    echo ""
    
    if [ $analysis_result -eq 0 ]; then
        echo -e "${GREEN}状态: ✓ 成功${NC}"
        return 0
    else
        echo -e "${RED}状态: ✗ 有部分失败${NC}"
        return 1
    fi
}

################################################################################
# 执行主程序
################################################################################

main "$@"
exit $?
