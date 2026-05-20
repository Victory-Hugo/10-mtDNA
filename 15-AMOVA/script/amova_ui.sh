#!/usr/bin/env bash
# ─────────────────────────────────────────────────────────────────────────────
# script/amova_ui.sh  —  AMOVA Pipeline 终端 UI 工具库
#
# 用法：source script/amova_ui.sh
# 所有 UI 函数写入 stderr，不影响 stdout 数据流。
# 在日志重定向（exec 2>&1）之前调用 amova_ui_init 可保证终端色彩检测正确。
# ─────────────────────────────────────────────────────────────────────────────

if [[ -n "${_AMOVA_UI_SH_LOADED:-}" ]]; then return 0; fi
_AMOVA_UI_SH_LOADED=1

# ══════════════════════════════════════════════════════════════════════════════
# 终端能力检测（幂等，多次调用安全）
# ══════════════════════════════════════════════════════════════════════════════

_AMOVA_UI_READY=0

amova_ui_init() {
    [[ "$_AMOVA_UI_READY" -eq 1 ]] && return 0
    _AMOVA_UI_READY=1

    local stderr_tty=0 color_count=0
    local locale="${LC_ALL:-${LC_CTYPE:-${LANG:-}}}"

    [[ -t 2 ]] && stderr_tty=1

    # ── 颜色支持 ──────────────────────────────────────────────────────────────
    _UI_COLOR=0
    if [[ -z "${NO_COLOR:-}" && "${TERM:-}" != "dumb" && "$stderr_tty" -eq 1 ]]; then
        if command -v tput >/dev/null 2>&1; then
            color_count=$(tput colors 2>/dev/null || printf '0')
        else
            color_count=8
        fi
        [[ "$color_count" =~ ^[0-9]+$ ]] && (( color_count >= 8 )) && _UI_COLOR=1
    fi

    # ── Unicode 支持 ───────────────────────────────────────────────────────────
    _UI_UNICODE=0
    case "$locale" in *UTF-8*|*utf8*|*utf-8*) _UI_UNICODE=1 ;; esac
    [[ -n "${AMOVA_ASCII_ONLY:-}" ]] && _UI_UNICODE=0

    # ── 色彩变量（青绿 / 金色 / 紫罗兰 配色方案）─────────────────────────────
    if [[ "$_UI_COLOR" -eq 1 ]]; then
        _C0=$'\033[0m'              # reset
        _CB=$'\033[1m'             # bold
        _CD=$'\033[2m'             # dim
        _CTEAL=$'\033[38;5;38m'    # 青绿（主色）
        _CGOLD=$'\033[38;5;214m'   # 金色（高亮数值）
        _CGRN=$'\033[38;5;119m'    # 嫩绿（OK / within）
        _CAMB=$'\033[38;5;208m'    # 琥珀（警告）
        _CRED=$'\033[38;5;196m'    # 红（错误）
        _CVIO=$'\033[38;5;141m'    # 紫罗兰（步骤编号 / among groups）
        _CSLV=$'\033[38;5;250m'    # 银灰（次要文字）
    else
        _C0="" _CB="" _CD="" _CTEAL="" _CGOLD=""
        _CGRN="" _CAMB="" _CRED="" _CVIO="" _CSLV=""
    fi

    # ── 图标集 ─────────────────────────────────────────────────────────────────
    if [[ "$_UI_UNICODE" -eq 1 ]]; then
        _IC_OK="✓"
        _IC_CACHED="⊙"
        _IC_WARN="△"
        _IC_ERR="✗"
        _IC_INFO="·"
        _IC_STEP="◆"
        _IC_ARROW="→"
        _IC_RULE="╌"
        _IC_BAR_FULL="█"
        _IC_BAR_LITE="░"
        _IC_TEE="├"
        _IC_END="└"
        _IC_VRT="│"
    else
        _IC_OK="+"
        _IC_CACHED="~"
        _IC_WARN="!"
        _IC_ERR="x"
        _IC_INFO="."
        _IC_STEP="#"
        _IC_ARROW="->"
        _IC_RULE="-"
        _IC_BAR_FULL="#"
        _IC_BAR_LITE="."
        _IC_TEE="|"
        _IC_END="'"
        _IC_VRT="|"
    fi
}

# ══════════════════════════════════════════════════════════════════════════════
# 内部工具
# ══════════════════════════════════════════════════════════════════════════════

_ui_rep() {
    # 重复字符 $1 共 $2 次
    local ch="$1" n="$2" out="" i=0
    while (( i < n )); do out+="$ch"; (( i++ )); done
    printf '%s' "$out"
}

_ui_hr() {
    # 水平分隔线（宽度默认 68）
    local w="${1:-68}"
    printf '%b%s%b\n' "$_CD" "$(_ui_rep "$_IC_RULE" "$w")" "$_C0" >&2
}

# ══════════════════════════════════════════════════════════════════════════════
# 消息级别
# ══════════════════════════════════════════════════════════════════════════════

ui_info() {
    amova_ui_init
    printf '%b%s%b  %s\n' "$_CTEAL" "$_IC_INFO" "$_C0" "$1" >&2
}

ui_ok() {
    amova_ui_init
    printf '%b%s%b  %s\n' "$_CGRN$_CB" "$_IC_OK" "$_C0" "$1" >&2
}

ui_cache_hit() {
    # Step 缓存命中（跳过重算）
    amova_ui_init
    local step="$1" title="$2" reason="${3:-指纹匹配}"
    printf '%b%s%b  Step %s cached  %b%s%b  %b%s%b\n' \
        "$_CTEAL$_CB" "$_IC_CACHED" "$_C0" \
        "$step" \
        "$_CSLV" "$title" "$_C0" \
        "$_CD" "$reason" "$_C0" >&2
}

ui_cache_miss() {
    # 缓存未命中或已失效
    amova_ui_init
    local reason="${1:-}"
    [[ -n "$reason" ]] && printf '%b%s%b  %s\n' "$_CAMB" "$_IC_INFO" "$_C0" "$reason" >&2
}

ui_warn() {
    amova_ui_init
    printf '%b%s%b  %s\n' "$_CAMB$_CB" "$_IC_WARN" "$_C0" "$1" >&2
}

ui_error() {
    amova_ui_init
    printf '%b%s%b  %s\n' "$_CRED$_CB" "$_IC_ERR" "$_C0" "$1" >&2
}

# ══════════════════════════════════════════════════════════════════════════════
# 步骤标题 / 完成
# ══════════════════════════════════════════════════════════════════════════════

ui_step_begin() {
    amova_ui_init
    local num="$1" title="$2"
    printf '\n' >&2
    _ui_hr
    printf '%b%s  Step %-2s%b  %b%s%b\n' \
        "$_CVIO$_CB" "$_IC_STEP" "$num" "$_C0" \
        "$_CB" "$title" "$_C0" >&2
    _ui_hr
}

ui_step_done() {
    amova_ui_init
    local num="$1" title="$2" secs="${3:-}"
    local time_hint=""
    [[ -n "$secs" ]] && time_hint="  ${_CD}elapsed ${secs}s${_C0}"
    printf '%b%s%b  Step %s done  %b%s%b%s\n' \
        "$_CGRN$_CB" "$_IC_OK" "$_C0" \
        "$num" \
        "$_CSLV" "$title" "$_C0" \
        "$time_hint" >&2
}

# ══════════════════════════════════════════════════════════════════════════════
# 键值 / 统计输出
# ══════════════════════════════════════════════════════════════════════════════

ui_kv() {
    # 通用键值对（配置参数展示）
    amova_ui_init
    local key="$1" val="$2"
    printf '  %b%-22s%b %b%s%b\n' \
        "$_CSLV" "$key" "$_C0" \
        "$_CB" "$val" "$_C0" >&2
}

ui_stat() {
    # AMOVA 统计结果（名称 + 金色高亮值）
    amova_ui_init
    local label="$1" value="$2"
    printf '  %b%s%b  %-28s%b%s%b\n' \
        "$_CTEAL" "$_IC_ARROW" "$_C0" \
        "$label" \
        "$_CGOLD$_CB" "$value" "$_C0" >&2
}

ui_check() {
    # 环境检查条目（pass / fail）
    amova_ui_init
    local item="$1" status="${2:-ok}"
    if [[ "$status" == "ok" ]]; then
        printf '  %b%s%b  %s\n' "$_CGRN$_CB" "$_IC_OK" "$_C0" "$item" >&2
    else
        printf '  %b%s%b  %s\n' "$_CRED$_CB" "$_IC_ERR" "$_C0" "$item" >&2
    fi
}

# ══════════════════════════════════════════════════════════════════════════════
# Banner / Logo
# ══════════════════════════════════════════════════════════════════════════════

ui_logo() {
    # 方差分量条形图：Va / Vb / Vc 三层，长度体现典型 AMOVA 权重关系
    amova_ui_init
    if [[ "$_UI_UNICODE" -eq 1 ]]; then
        printf '  %b%-18s%b %b%s%b%s\n' \
            "$_CSLV" "Va  among groups" "$_C0" \
            "$_CVIO$_CB" "$(_ui_rep "$_IC_BAR_FULL" 5)" "$_C0" \
            "${_CD}$(_ui_rep "$_IC_BAR_LITE" 19)${_C0}" >&2
        printf '  %b%-18s%b %b%s%b%s\n' \
            "$_CSLV" "Vb  among pops" "$_C0" \
            "$_CTEAL$_CB" "$(_ui_rep "$_IC_BAR_FULL" 10)" "$_C0" \
            "${_CD}$(_ui_rep "$_IC_BAR_LITE" 14)${_C0}" >&2
        printf '  %b%-18s%b %b%s%b%s\n' \
            "$_CSLV" "Vc  within pops" "$_C0" \
            "$_CGRN$_CB" "$(_ui_rep "$_IC_BAR_FULL" 24)" "$_C0" \
            "" >&2
        printf '  %b%s%b\n' "$_CD" \
            "$(_ui_rep ' ' 20)└── FST = (Va+Vb) / (Va+Vb+Vc)" "$_C0" >&2
    else
        printf '  %-18s %s\n' "Va  among groups"  "##### ..................." >&2
        printf '  %-18s %s\n' "Vb  among pops"    "########## .............." >&2
        printf '  %-18s %s\n' "Vc  within pops"   "########################" >&2
        printf '  %s\n' "                    FST = (Va+Vb)/(Va+Vb+Vc)" >&2
    fi
}

ui_banner() {
    # 启动 banner：logo + 标题 + 版本信息
    amova_ui_init
    printf '\n' >&2
    _ui_hr
    ui_logo
    printf '\n' >&2
    printf '  %bAMOVA Pipeline%b\n' "$_CB" "$_C0" >&2
    printf '  %bAnalysis of Molecular Variance%b\n' "$_CTEAL" "$_C0" >&2
    printf '  %bExcoffier et al. (1992) Genetics 131:479-491%b\n' "$_CD" "$_C0" >&2
    _ui_hr
    printf '\n' >&2
}

ui_final_summary() {
    # 流程结束时的汇总框
    amova_ui_init
    printf '\n' >&2
    _ui_hr
    printf '  %b%s%b  %bPipeline complete%b\n' \
        "$_CGRN$_CB" "$_IC_OK" "$_C0" "$_CB" "$_C0" >&2
    _ui_hr
}
