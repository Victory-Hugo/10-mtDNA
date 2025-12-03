#!/bin/bash
set -euo pipefail

############################################
# 路径配置
############################################
#! 脚本使用R代码，需要提前安装 install.packages("fitdistrplus")
LIST_FILE="/mnt/c/Users/Administrator/Desktop/example/conf/Sample_list.tsv"     # 样本清单（两列：BAM路径\tSampleID）
OUTPUT_DIR="/mnt/c/Users/Administrator/Desktop/example/output"             # 输出 fasta 存放位置
TEMP_DIR="/mnt/c/Users/Administrator/Desktop/example/temp"                 # 中间文件存放位置
LOG_FILE="/mnt/c/Users/Administrator/Desktop/example/log/1_schmutzi_aDNA.log"

############################################
# 运行参数配置
############################################
export THREADS=8                # bwa / schmutzi 使用线程数
export LENGTH_DEAM=2            # contDeam length (如果失败可尝试改为 2 或 3)
export LIB_TYPE="double"        # contDeam library (单端文库用 "single")
export MIN_READS=50             # 跳过 deamination 的最小 reads 阈值
export USE_ENDO=0               # 是否使用 -endo 选项 (0=关闭, 1=开启)
                                # 如果样本损伤较低或出现 NaN 错误，设为 0
export ITERATIONS=2             # schmutzi 迭代次数 (1-5, 默认 2)
export LOG2FASTA_QUAL=20        # log2fasta 的质量阈值 (-q)
export LOG2FASTA_INDEL=0        # log2fasta 的 indel 参数 (-indel)
export RESUME=1                 # 断点续跑 (0=从头开始, 1=跳过已完成的样本)

############################################
# 软件路径配置
############################################
SOFTWARE_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/8-schmutzi"
export REF="$SOFTWARE_DIR/share/schmutzi/refs/human_MT.fa"
export AL_FREQ_DIR="$SOFTWARE_DIR/share/schmutzi/alleleFreqMT/eurasian/freqs"
export SCHMUTZI_PERL="$SOFTWARE_DIR/src/schmutzi.pl"
export LOG2FASTA="$SOFTWARE_DIR/src/log2fasta"
SCHMUTZI="${SOFTWARE_DIR}/script/1_schmutzi_aDNA.sh"
export SAMTOOLS_PATH="samtools"    # 如果samtools在PATH中，直接写samtools；否则写完整路径
export BWA_PATH="bwa"              # 如果bwa在PATH中，直接写bwa；否则写完整路径

mkdir -p "$OUTPUT_DIR" "$TEMP_DIR"

"${SCHMUTZI}" \
  "$LIST_FILE" \
  "$OUTPUT_DIR" \
  "$TEMP_DIR" \
  "$SOFTWARE_DIR" 

echo "任务已后台启动，日志：$LOG_FILE"
