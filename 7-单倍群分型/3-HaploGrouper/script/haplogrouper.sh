#!/usr/bin/bash

# 设置脚本的错误处理
set -euo pipefail

# 脚本主目录配置
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$(dirname "$(dirname "${SCRIPT_DIR}")")")"

# ============ 本地软件配置 ============
PYTHON="python3"
EXTERNAL_DATA="/mnt/e/Scientifc_software/haploGrouper"

# ============ 项目路径配置 ============
SRC_DIR="${PROJECT_ROOT}/3-HaploGrouper/src"
HGR="${SRC_DIR}/hGrpr2.py"
DATA_DIR="${PROJECT_ROOT}/3-HaploGrouper/data"
INPUT_DIR="${PROJECT_ROOT}/3-HaploGrouper/input"

# ============ 输入文件配置 ============
INPUT_VCF="${INPUT_DIR}/AltaicGuizhou_chr26.vcf"
TREE_FILE="${EXTERNAL_DATA}/data/mt_phyloTree_b17_Tree2.txt"
MUTATION_FILE="${EXTERNAL_DATA}/data/mt_phyloTree_b17_Mutation.txt"
REF_FASTA="${EXTERNAL_DATA}/data/rCRS.fasta"

# ============ 输出文件配置 ============
PREFIX="AltaicGuizhou_chr26.id_hGrpr2"
OUTPUT_HG="AltaicGuizhou_chr26_mt_hg_hGrpr2.txt"
OUTPUT_SCORES="AltaicGuizhou_chr26_mt_allScores_hGrpr2.txt"

# ============ 验证必要文件是否存在 ============
for file in "${HGR}" "${INPUT_VCF}" "${TREE_FILE}" "${MUTATION_FILE}" "${REF_FASTA}"; do
    if [[ ! -f "$file" ]]; then
        echo "错误：必要文件不存在 - $file" >&2
        exit 1
    fi
done

echo "开始运行HaploGrouper分析..."
echo "输入文件：${INPUT_VCF}"

# ============ 运行hGrpr2.py ============
# 使用变量替代硬编码路径，确保参数与参数之间用空格隔开
${PYTHON} "${HGR}" \
    -v "${INPUT_VCF}" \
    -t "${TREE_FILE}" \
    -l "${MUTATION_FILE}" \
    -f "${REF_FASTA}" \
    -i "${PREFIX}" \
    -o "${OUTPUT_HG}" \
    -x "${OUTPUT_SCORES}"

echo "分析完成！"
echo "输出文件："
echo "  - ${OUTPUT_HG}"
echo "  - ${OUTPUT_SCORES}"  