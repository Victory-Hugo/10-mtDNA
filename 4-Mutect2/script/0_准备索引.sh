#!/bin/bash
set -euo pipefail  # 严格模式：任何错误都终止脚本

# ------------------ 读取配置文件 ------------------
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONF_FILE="${SCRIPT_DIR}/../conf/2_Mutect2.conf"
if [[ ! -f "${CONF_FILE}" ]]; then
    echo "Config not found: ${CONF_FILE}" >&2
    exit 1
fi
source "${CONF_FILE}"

# ------------------ 文件检查与清理 ------------------
echo "清理旧的索引文件和字典文件（如存在）..."

[ -f "${REF_CHRM}.fai" ] && rm -f "${REF_CHRM}.fai"
[ -f "${REF_CHRM_SHIFTED}.fai" ] && rm -f "${REF_CHRM_SHIFTED}.fai"

[ -f "${REF_CHRM_DICT}" ] && rm -f "${REF_CHRM_DICT}"
[ -f "${REF_CHRM_SHIFTED_DICT}" ] && rm -f "${REF_CHRM_SHIFTED_DICT}"

[ -f "${REF_CHRM_SHIFTED}" ] && rm -f "${REF_CHRM_SHIFTED}"
[ -f "${REF_SHIFT_CHAIN}" ] && rm -f "${REF_SHIFT_CHAIN}"

# ------------------ 创建原始参考索引 ------------------
echo "创建原始参考序列索引..."
${SAMTOOLS} faidx "${REF_CHRM}"

echo "创建原始参考序列字典..."
${PICARD} CreateSequenceDictionary \
    R="${REF_CHRM}" \
    O="${REF_CHRM_DICT}"

# ------------------ GATK 生成 shifted FASTA ------------------
echo "生成 shifted FASTA 及 shift back chain 文件..."
${GATK} ShiftFasta \
  -R "${REF_CHRM}" \
  -O "${REF_CHRM_SHIFTED}" \
  --shift-back-output "${REF_SHIFT_CHAIN}"

# ------------------ 创建 shifted 参考索引 ------------------
echo "创建 shifted FASTA 的 BWA 索引..."
${BWA} index "${REF_CHRM_SHIFTED}"

echo "创建 shifted FASTA 的字典..."
${PICARD} CreateSequenceDictionary \
    R="${REF_CHRM_SHIFTED}" \
    O="${REF_CHRM_SHIFTED_DICT}"

echo "所有参考索引和字典准备完成。"
