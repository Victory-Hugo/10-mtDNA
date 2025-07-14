#!/bin/bash
set -euo pipefail  # 严格模式：任何错误都终止脚本

# ------------------ 路径变量定义 ------------------
REF_DIR="/mnt/e/Scientifc_software/TT/Reference"
GATK="/mnt/e/Scientifc_software/gatk-4.4.0.0/gatk"

ORIG_FASTA="${REF_DIR}/chrM_rCRS.fasta"
SHIFTED_FASTA="${REF_DIR}/chrM_rCRS.shifted.fa"
SHIFT_CHAIN="${REF_DIR}/chrM_rCRS.shift_back.chain"

ORIG_DICT="${REF_DIR}/chrM_rCRS.dict"
SHIFTED_DICT="${REF_DIR}/chrM_rCRS.shifted.dict"

# ------------------ 文件检查与清理 ------------------
echo "清理旧的索引文件和字典文件（如存在）..."

[ -f "${ORIG_FASTA}.fai" ] && rm -f "${ORIG_FASTA}.fai"
[ -f "${SHIFTED_FASTA}.fai" ] && rm -f "${SHIFTED_FASTA}.fai"

[ -f "${ORIG_DICT}" ] && rm -f "${ORIG_DICT}"
[ -f "${SHIFTED_DICT}" ] && rm -f "${SHIFTED_DICT}"

[ -f "${SHIFTED_FASTA}" ] && rm -f "${SHIFTED_FASTA}"
[ -f "${SHIFT_CHAIN}" ] && rm -f "${SHIFT_CHAIN}"

# ------------------ 创建原始参考索引 ------------------
echo "创建原始参考序列索引..."
samtools faidx "${ORIG_FASTA}"

echo "创建原始参考序列字典..."
picard CreateSequenceDictionary \
    R="${ORIG_FASTA}" \
    O="${ORIG_DICT}"

# ------------------ GATK 生成 shifted FASTA ------------------
echo "生成 shifted FASTA 及 shift back chain 文件..."
"${GATK}" ShiftFasta \
  -R "${ORIG_FASTA}" \
  -O "${SHIFTED_FASTA}" \
  --shift-back-output "${SHIFT_CHAIN}"

# ------------------ 创建 shifted 参考索引 ------------------
echo "创建 shifted FASTA 的 BWA 索引..."
bwa index "${SHIFTED_FASTA}"

echo "创建 shifted FASTA 的字典..."
picard CreateSequenceDictionary \
    R="${SHIFTED_FASTA}" \
    O="${SHIFTED_DICT}"

echo "✅ 所有参考索引和字典准备完成。"
