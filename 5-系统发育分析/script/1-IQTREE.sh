#!/bin/bash

IQTREE3_BIN="/mnt/e/Scientifc_software/iqtree-3.0.0-Linux-intel/bin/iqtree3"
INPUT_FASTA="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/5-系统发育分析/example/LLST_500.aln.fasta"
OUTPUT_PREFIX="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/5-系统发育分析/output/LLST_500"

 ${IQTREE3_BIN} \
    -s ${INPUT_FASTA} \
    -m "MIX{GTR+FO,GTR+FO,GTR+FO,HKY+FO}+I+R5" \
    -pre ${OUTPUT_PREFIX} \
    -T AUTO

# halign4 \
#     "/mnt/c/Users/Administrator/Desktop/LLST_500.fasta" \
#     "/mnt/c/Users/Administrator/Desktop/LLST_500.aln.fasta" \
#     -r "/mnt/f/OneDrive/文档（科研）/脚本/Download/2-Dimensionality-analysis/4-DAPC/conf/chrM.fasta" \
#     -t 16 \
#     -sa 15


