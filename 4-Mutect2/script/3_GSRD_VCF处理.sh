#!/usr/bin/env bash
set -euo pipefail

# 配置目录
VCF_DIR="/mnt/l/22-WHALE/1-数据/data/低深度/Final" # 输入目录，包含待处理的 .vcf / .vcf.gz 文件
TEMP_DIR="/mnt/l/22-WHALE/1-数据/data/低深度/temp" # 临时文件目录,随便写
OUT_DIR="/mnt/l/22-WHALE/1-数据/data/低深度/output" # 输出目录
PYTHON="/home/luolintao/miniconda3/envs/BigLin/bin/python3" # Python 解释器路径
SCRIPT="/mnt/f/onedrive/文档（科研）/脚本/Download/10-mtDNA/4-Mutect2/script/3_GSRD_VCF处理.py"

# 创建输出／临时目录
mkdir -p "$TEMP_DIR" "$OUT_DIR"

export PYTHON SCRIPT TEMP_DIR OUT_DIR

# 并行处理每个 .vcf / .vcf.gz
find "$VCF_DIR" -maxdepth 1 -type f \( -name '*.vcf' -o -name '*.vcf.gz' \) | \
  parallel --bar -j 8 '
    f={}
    base=$(basename "$f")
    base=${base%.vcf.gz}
    base=${base%.vcf}
    tmp="$TEMP_DIR/${base}.filtered_clean.vcf"
    out="$OUT_DIR/${base}.filtered_clean.vcf"

    # 1) awk 清洗，兼容未压缩/压缩 VCF 输入
    if [[ "$f" == *.vcf.gz ]]; then
      gzip -dc -- "$f"
    else
      cat -- "$f"
    fi | awk '\''BEGIN{FS=OFS="\t"}
      /^#/ { print; next }
      $7=="PASS" {
        key=$1":"$2
        if (!(key in seen)) {
          seen[key]=1
          $8="."
          print
        }
      }
    '\'' > "$tmp"

    # 2) python 脚本进一步处理
    $PYTHON $SCRIPT -i "$tmp" -o "$out"
  '

# 清除临时文件
rm -rf "$TEMP_DIR"

echo "All done. Cleaned up $TEMP_DIR."
