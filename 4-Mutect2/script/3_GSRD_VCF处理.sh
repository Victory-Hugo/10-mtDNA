#!/usr/bin/env bash
set -euo pipefail

# 配置目录
VCF_DIR="/mnt/e/Scientifc_software/Haplogrep3/input"
TEMP_DIR="/mnt/e/Scientifc_software/Haplogrep3/temp"
OUT_DIR="/mnt/e/Scientifc_software/Haplogrep3/input"
PYTHON="/home/luolintao/miniconda3/bin/python3"
SCRIPT="/mnt/e/Scientifc_software/TT/script/3_GSRD_VCF处理.py"

# 创建输出／临时目录
mkdir -p "$TEMP_DIR" "$OUT_DIR"

export PYTHON SCRIPT TEMP_DIR OUT_DIR

# 并行处理每个 .vcf
find "$VCF_DIR" -maxdepth 1 -type f -name '*.vcf' | \
  parallel --bar -j 8 '
    f={}
    base=$(basename "$f" .vcf)
    tmp="$TEMP_DIR/${base}.filtered_clean.vcf"
    out="$OUT_DIR/${base}.filtered_clean.vcf"

    # 1) awk 清洗
    awk '\''BEGIN{FS=OFS="\t"}
      /^#/ { print; next }
      $7=="PASS" {
        key=$1":"$2
        if (!(key in seen)) {
          seen[key]=1
          $8="."
          print
        }
      }
    '\'' "$f" > "$tmp"

    # 2) python 脚本进一步处理
    $PYTHON $SCRIPT -i "$tmp" -o "$out"
  '

# 清除临时文件
rm -rf "$TEMP_DIR"

echo "All done. Cleaned up $TEMP_DIR."
