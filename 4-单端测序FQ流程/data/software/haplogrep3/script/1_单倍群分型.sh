#!/bin/bash

#*###############################
#*########### VCF #############
#*###############################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOFTWARE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

# 设置输入和输出目录
input_dir="${SOFTWARE_DIR}/examples/input/vcf"
output_dir="${SOFTWARE_DIR}/examples/output"

# 检查输出目录是否存在，如果不存在则创建
mkdir -p "$output_dir"
cd ${SOFTWARE_DIR}
# 遍历目录中的所有 vcf 文件
for vcf_file in "$input_dir"/*.vcf; do
    # 获取不带路径的文件名
    filename=$(basename "$vcf_file")
    # 构建输出文件名
    output_file="${output_dir}/${filename%.vcf}.txt"

    # 执行 classify 命令
    ${SOFTWARE_DIR}/haplogrep3 classify \
        --input "$vcf_file" \
        --tree phylotree-rcrs@17.2 \
        --output "$output_file" \
        --extend-report \
        --write-fasta  #! 如果是芯片加上--chip

    # 输出处理的文件名
    echo "Processed: $vcf_file"
done

# #*###############################
# #*########### FASTA #############
# #*###############################

# # 设置输入和输出目录
# input_dir="${SOFTWARE_DIR}/examples/input/fasta"
# output_dir="${SOFTWARE_DIR}/examples/output"

# # 检查输出目录是否存在，如果不存在则创建
# mkdir -p "$output_dir"
# cd ${SOFTWARE_DIR}
# # 遍历目录中的所有 fasta 文件
# for fasta_file in "$input_dir"/*.fasta; do
#     # 获取不带路径的文件名
#     filename=$(basename "$fasta_file")
#     # 构建输出文件名
#     output_file="${output_dir}/${filename%.vcf}.txt"

#     # 执行 classify 命令
#     ${SOFTWARE_DIR}/haplogrep3 classify \
#         --input "$fasta_file" \
#         --tree phylotree-rcrs@17.2 \
#         --output "$output_file" \
#         --extend-report 
#     # 输出处理的文件名
#     echo "Processed: $fasta_file"
# done


# #*###############################
# #*########### 单倍群聚类 #########
# #*###############################

# ./haplogrep3 cluster-haplogroups \
# 	--output "${SOFTWARE_DIR}/examples/output" \
# 	--tree phylotree-rcrs@17.2


# #*###############################
# #*########### 序列比对 #########
# #*###############################


# ./haplogrep3 aligned \
# 	--fasta /All.fasta \
# 	--output /HGDP.aln.fasta \
# 	--tree phylotree-rcrs@17.2
