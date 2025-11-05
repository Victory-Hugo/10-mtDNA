#!/bin/bash

DATA_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/7-单倍群分型/1-HaploCart/data" #? 数据存放位置
cd $DATA_DIR

# 原版数据下载
wget -nc -l1 \
    --recursive --no-directories --no-parent \
    -P $DATA_DIR \
    ftp://ftp.healthtech.dtu.dk:/public/haplocart/hcfiles/
wget ftp://ftp.healthtech.dtu.dk:/public/haplocart/testdata/rCRS.fa.gz

#?↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑↑
# 123网盘转存2025年11月5日留档
# https://www.123865.com/s/9k71Td-W5sJ3

#* ======基础命令============
# # 对于共识 FASTA（可以压缩）：
# vgan haplocart -f rCRS.fa.gz

# # 对于单端交错（可以压缩）：
# wget ftp://ftp.healthtech.dtu.dk:/public/haplocart/testdata/rCRS.fq.gz
# vgan haplocart -fq1 rCRS.fq.gz

# # 对于配对端 FASTQ（可以压缩）：
# vgan haplocart -fq1 [FASTQ file] -fq2 [FASTQ file]