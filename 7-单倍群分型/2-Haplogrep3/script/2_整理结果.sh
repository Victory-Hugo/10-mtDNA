#!/usr/bin/env bash
set -euo pipefail

# 目标目录
OUTPUT_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/7-单倍群分型/2-Haplogrep3/examples/output"

# 最终合并文件
QC_ALL="$OUTPUT_DIR/All.qc.txt"
TXT_ALL="$OUTPUT_DIR/All.txt"
FASTA_ALL="$OUTPUT_DIR/All.fasta"

# 清空旧的结果文件（若存在）
: > "$QC_ALL"
: > "$TXT_ALL"
: > "$FASTA_ALL"

# 1) 合并 .qc.txt，只保留一个标题行
first=1
for f in "$OUTPUT_DIR"/*.qc.txt; do
  [ ! -s "$f" ] && continue
  if [ $first -eq 1 ]; then
    cat "$f" >> "$QC_ALL"
    first=0
  else
    tail -n +2 "$f" >> "$QC_ALL"
  fi
done

# 2) 合并 .txt（排除 .qc.txt），只保留一个标题行
first=1
for f in "$OUTPUT_DIR"/*.txt; do
  [[ "$f" == *.qc.txt ]] && continue
  [ ! -s "$f" ] && continue
  if [ $first -eq 1 ]; then
    cat "$f" >> "$TXT_ALL"
    first=0
  else
    tail -n +2 "$f" >> "$TXT_ALL"
  fi
done

# 3) 合并所有 .fasta
for f in "$OUTPUT_DIR"/*.fasta; do
  cat "$f" >> "$FASTA_ALL"
done

# 4) 删除原始单文件，只保留 All.qc.txt、All.txt、All.fasta
for f in "$OUTPUT_DIR"/*.{qc.txt,txt,fasta}; do
  # 跳过合并结果文件
  if [[ "$f" != "$QC_ALL" && "$f" != "$TXT_ALL" && "$f" != "$FASTA_ALL" ]]; then
    rm -f "$f"
  fi
done

echo "Done. 只保留合并文件："
echo "  $QC_ALL"
echo "  $TXT_ALL"
echo "  $FASTA_ALL"
