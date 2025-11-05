#!/usr/bin/env bash
# ============================================================
# Author: BigLin
# Purpose: 将多序列 FASTA 拆分为单序列 FASTA；文件名=序列ID（不含描述）
# Note: 变量全部写死。包含彩色打印与进度条。
# ============================================================

# -------- 写死的配置（按需自行修改后再运行） -------------------
INPUT="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/7-单倍群分型/1-HaploCart/example/example.fasta"      # 多序列 FASTA 输入文件
OUTDIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/7-单倍群分型/1-HaploCart/example/"     # 输出目录
BAR_WIDTH=40                                   # 进度条宽度（字符数）
# ------------------------------------------------------------

# 颜色（ANSI 转义）
C_RESET="\033[0m"
C_RED="\033[31m"
C_GREEN="\033[32m"
C_YELLOW="\033[33m"
C_BLUE="\033[34m"
C_CYAN="\033[36m"

# 打印函数
info()  { printf "${C_BLUE}[INFO]${C_RESET} %s\n" "$*"; }
warn()  { printf "${C_YELLOW}[WARN]${C_RESET} %s\n" "$*"; }
ok()    { printf "${C_GREEN}[OK]${C_RESET}   %s\n" "$*"; }
fail()  { printf "${C_RED}[FAIL]${C_RESET} %s\n" "$*"; }

# 进度条函数
# 参数: 已完成数 总数
draw_progress() {
  local done=$1 total=$2
  local percent=0
  if [ "$total" -gt 0 ]; then
    percent=$(( done * 100 / total ))
  fi
  local filled=$(( done * BAR_WIDTH / total ))
  local empty=$(( BAR_WIDTH - filled ))
  printf "\r${C_CYAN}["
  printf "%0.s#" $(seq 1 $filled)
  printf "%0.s-" $(seq 1 $empty)
  printf "] %3d%% (${done}/${total})${C_RESET}" "$percent"
}

# 预检
if [ ! -f "$INPUT" ]; then
  fail "找不到输入文件: $INPUT"
  exit 1
fi

mkdir -p "$OUTDIR" || { fail "无法创建输出目录: $OUTDIR"; exit 1; }

# 统计序列数量（以 '>' 开头的行）
TOTAL=$(grep -c '^>' "$INPUT" 2>/dev/null || echo 0)
if [ "$TOTAL" -eq 0 ]; then
  fail "未在输入文件中检测到任何序列头（以 '>' 开头）: $INPUT"
  exit 1
fi

info "输入文件: $INPUT"
info "输出目录: $OUTDIR"
info "检测到序列数: $TOTAL"
info "开始拆分..."

# 临时管道：awk 输出每个完成的序列 ID，bash 侧更新进度条
# 说明：
# - 仅使用头部的第一个字段（遇到空白停止）作为文件名
# - 原始序列按输入保留；输出头仅为 >ID（不含描述）
# - 对文件名做轻度清洗（替换 / \ : * ? \" < > | 和空白 为 下划线）
PROCESSED=0
TMPLOG="$(mktemp)"
trap 'rm -f "$TMPLOG"' EXIT

# 用 awk 拆分写文件，同时将 ID 打印到 stdout 供进度条使用
awk -v outdir="$OUTDIR" '
  BEGIN{
    RS=">"; ORS="";
  }
  NR>1 {
    # 取 header（第一行）与序列主体
    header = $1
    sub(/\r/,"", header)
    # 提取ID（至第一个空白）
    split(header, a, /[ \t\r\n]+/)
    id = a[1]

    # 清洗为安全文件名（与 bash 保持一致的字符集合）
    gsub(/[\/\\:\*\?"<>\|\t ]/, "_", id)

    # 获取序列主体（去掉 header 行）
    body = substr($0, index($0, "\n") + 1)
    gsub(/\r/, "", body)

    # 输出文件
    fname = outdir "/" id ".fasta"
    print ">" id "\n" body > fname
    close(fname)

    # 将完成的 id 打印给外层进度条
    print id "\n" > "/dev/fd/3"
  }
' 3> "$TMPLOG" "$INPUT" &

# 读取进度并绘制条形
# 说明：通过 tail -f + 等待 awk 持续写入，再逐行更新
# 为了简洁稳健，这里使用 while 循环读取新行
# 注：不使用 set -e 以避免管道早退
{
  # 循环直到 awk 结束且文件读尽
  while true; do
    # 逐行读取新生成的 ID
    while IFS= read -r _line; do
      PROCESSED=$((PROCESSED + 1))
      draw_progress "$PROCESSED" "$TOTAL"
    done < "$TMPLOG"

    # 若后台 awk 仍在运行，短暂等待并重试；否则跳出
    if ps -p $! >/dev/null 2>&1; then
      sleep 0.1
      continue
    else
      # 再扫尾一次，避免竞态遗漏
      while IFS= read -r _line; do
        PROCESSED=$((PROCESSED + 1))
        draw_progress "$PROCESSED" "$TOTAL"
      done < "$TMPLOG"
      break
    fi
  done
} >/dev/tty 2>/dev/null || true

# 换行收尾
echo

# 结束提示
if [ "$PROCESSED" -eq "$TOTAL" ]; then
  ok "拆分完成。共生成 ${PROCESSED} 个 FASTA 文件。"
  info "示例文件: $(ls -1 "$OUTDIR" | head -n 3 | xargs -I{} echo "$OUTDIR/{}")"
else
  warn "拆分未完整结束：期望 $TOTAL，实际 $PROCESSED。"
fi
