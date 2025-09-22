#!/bin/bash

FST_R_SCRIPT="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/5-多样性和Fst分析/Fst/src/0-Fst.R"
INPUT_ALN_FASAT="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/5-多样性和Fst分析/Fst/example/Example.aln.fasta"
OUTPUT_FST_MATRIX="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/5-多样性和Fst分析/Fst/output/fst_matrix.csv"
POPULATION_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/5-多样性和Fst分析/Fst/example/Example.txt"

# 可选的Fst计算方法：
# - hudson    : Hudson et al. (1992) F_ST (默认)
# - nei       : Nei's G_ST 
# - nucleotide: 基于核苷酸的F_ST
# - haplotype : 基于单倍型的F_ST
# Rscript \  
#     src/0-Fst.R \
#     example/Example.aln.fasta \
#     example/Example.txt \
#     output/fst_basic.csv \
#     hudson \ #* 选用的Fst计算方法
#     16 \ #* 使用的CPU核心数
#     0 \ #* bootstrap重采样次数（设为0跳过统计检验）
#     0.05 #* 显著性水平

# 只计算Fst值（无统计检验）
Rscript \
    src/0-Fst.R \
    example/Example.aln.fasta \
    example/Example.txt \
    output/fst_basic.csv \
    hudson \
    16 \
    0 \
    0.05

# 删除中间文件
rm -rf temp_popgenome*