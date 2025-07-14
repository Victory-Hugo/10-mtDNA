: '
脚本名称: 2_Mutect2.sh

用途说明:
本脚本用于批量处理线粒体DNA (mtDNA) 的变异检测流程，适用于多个样本的自动化分析。主要流程包括：
1. 参考基因组索引建立（如首次运行需取消注释）。
2. 读取包含样本名和BAM文件路径的列表文件，逐个样本并行处理。
3. 对每个样本，完成chrM提取、排序、索引、变异检测、重比对、去重、指标收集、变异整合与过滤等步骤。
4. 支持GNU parallel并行加速。

参数说明:
- 输入参数: 仅需一个参数，为包含样本名和BAM文件路径的列表文件（每行两列，制表符或空格分隔）。
- 依赖工具: samtools, picard, gatk, bwa, bcftools, bgzip, GNU parallel
- 参考文件: 需提前准备好参考基因组、字典、shifted参考、chain文件等。

使用方法:
bash 2_Mutect2.sh 4-Mutect2/meta/list.txt

输出说明:
每个样本在指定输出目录下生成一系列中间文件和最终的变异检测结果（VCF格式）。

注意事项:
- 首次运行前请确保参考基因组已建立索引和字典。
- 需根据实际环境调整路径和线程数。
- 需保证所有依赖工具已正确安装并在PATH中。
'
#!/bin/bash
set -euo pipefail

# ==================================================
# 基础路径定义
BASE_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/4-Mutect2"
REF_DIR="${BASE_DIR}/conf"
OUT_DIR="${BASE_DIR}/output"  # 基础输出目录

# 参考基因组及其字典定义
WHOLE_GENOME_REFERENCE="${REF_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
WHOLE_GENOME_REFERENCE_DICT="${REF_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta.dict"

# 通用参数定义
JAVA_OPT="-Xmx64g -DGATK_STACKTRACE_ON_USER_EXCEPTION=true"
THREADS=16
# ==================================================

# 如果需要建立参考基因组索引（只需执行一次），可以取消下面注释
samtools faidx "${REF_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
picard CreateSequenceDictionary \
    R="${WHOLE_GENOME_REFERENCE}" \
    O="${WHOLE_GENOME_REFERENCE_DICT}"

# 检查传入参数：必须传入含有两列的文件
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 list.txt"
    exit 1
fi

LIST_FILE="$1"

# 单个样本处理函数
process_sample() {
    local base_name bam_file
    base_name="$1"
    bam_file="$2"

    echo "Processing sample: ${base_name}"

    # 为当前样本创建输出目录
    SAMPLE_OUT="${OUT_DIR}/${base_name}"
    mkdir -p "${SAMPLE_OUT}"

    # ---- 原始数据处理 ----
    gatk --java-options "${JAVA_OPT}" PrintReads \
        -I "$bam_file" \
        -L chrM \
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter \
        -O "${SAMPLE_OUT}/${base_name}.chrM.bam"

    samtools sort -@ ${THREADS} -O bam \
        -o "${SAMPLE_OUT}/${base_name}.chrM.sorted.bam" \
        "${SAMPLE_OUT}/${base_name}.chrM.bam"

    picard BuildBamIndex \
        INPUT="${SAMPLE_OUT}/${base_name}.chrM.sorted.bam" \
        OUTPUT="${SAMPLE_OUT}/${base_name}.chrM.sorted.bai"

    # ---- 非控制区变异检测 ----
    gatk --java-options "${JAVA_OPT}" Mutect2 \
        -R "${REF_DIR}/chrM_rCRS.fasta" \
        -L chrM \
        --mitochondria-mode \
        --annotation StrandBiasBySample \
        --max-reads-per-alignment-start 75 \
        --max-mnp-distance 0 \
        -I "${SAMPLE_OUT}/${base_name}.chrM.sorted.bam" \
        -O "${SAMPLE_OUT}/${base_name}.chrM.ncr.vcf"

    sed '/^##contig.*\(random\|Un\|_alt\|HLA\|EBV\|chr[1-9XY]\)/d' \
        "${SAMPLE_OUT}/${base_name}.chrM.ncr.vcf" \
        > "${SAMPLE_OUT}/${base_name}.chrM.ncr.reform.vcf"

    bgzip -f -@ ${THREADS} "${SAMPLE_OUT}/${base_name}.chrM.ncr.reform.vcf"
    bcftools index --threads ${THREADS} -f -t \
        "${SAMPLE_OUT}/${base_name}.chrM.ncr.reform.vcf.gz"

    # ---- 数据重处理流程 ----
    gatk RevertSam \
        -I "${SAMPLE_OUT}/${base_name}.chrM.sorted.bam" \
        -O "${SAMPLE_OUT}/${base_name}.chrM.reverted.bam"

    picard SamToFastq \
        I="${SAMPLE_OUT}/${base_name}.chrM.reverted.bam" \
        FASTQ="${SAMPLE_OUT}/${base_name}_R1.chrM.fq.gz" \
        SECOND_END_FASTQ="${SAMPLE_OUT}/${base_name}_R2.chrM.fq.gz"

    # ---- 重比对流程 ----
    bwa mem -t ${THREADS} -M \
        -R "@RG\tID:group1\tSM:${base_name}\tPL:illumina\tLB:CasCADE\tPU:unit1" \
        "${REF_DIR}/chrM_rCRS.shifted.fa" \
        "${SAMPLE_OUT}/${base_name}_R1.chrM.fq.gz" \
        "${SAMPLE_OUT}/${base_name}_R2.chrM.fq.gz" \
        > "${SAMPLE_OUT}/${base_name}.chrM.shifted.sam"

    samtools sort -@ ${THREADS} -O bam \
        -o "${SAMPLE_OUT}/${base_name}.chrM.shifted.sorted.bam" \
        "${SAMPLE_OUT}/${base_name}.chrM.shifted.sam"

    # ---- 数据标记和指标收集 ----
    picard AddOrReplaceReadGroups \
        I="${SAMPLE_OUT}/${base_name}.chrM.shifted.sorted.bam" \
        O="${SAMPLE_OUT}/${base_name}.chrM.shifted.rg.bam" \
        RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM="${base_name}"

    picard MarkDuplicates \
        INPUT="${SAMPLE_OUT}/${base_name}.chrM.shifted.rg.bam" \
        OUTPUT="${SAMPLE_OUT}/${base_name}.chrM.dedup.bam" \
        METRICS_FILE="${SAMPLE_OUT}/${base_name}.chrM.dedup.metrics"

    picard BuildBamIndex \
        INPUT="${SAMPLE_OUT}/${base_name}.chrM.dedup.bam" \
        OUTPUT="${SAMPLE_OUT}/${base_name}.chrM.dedup.bai"

    picard CollectWgsMetrics \
        I="${SAMPLE_OUT}/${base_name}.chrM.dedup.bam" \
        O="${SAMPLE_OUT}/${base_name}.chrM.wgs_metrics.txt" \
        R="${REF_DIR}/chrM_rCRS.shifted.fa"

    # ---- 控制区变异检测 ----
    gatk Mutect2 \
        -R "${REF_DIR}/chrM_rCRS.shifted.fa" \
        -L chrM \
        --mitochondria-mode \
        --annotation StrandBiasBySample \
        --max-reads-per-alignment-start 75 \
        --max-mnp-distance 0 \
        -I "${SAMPLE_OUT}/${base_name}.chrM.dedup.bam" \
        -O "${SAMPLE_OUT}/${base_name}.chrM.cr.vcf.gz"

    # ---- 变异数据整合 ----
    picard LiftoverVcf \
        I="${SAMPLE_OUT}/${base_name}.chrM.cr.vcf.gz" \
        O="${SAMPLE_OUT}/${base_name}.chrM.cr.liftover.vcf.gz" \
        CHAIN="${REF_DIR}/chrM_rCRS.shift_back.chain" \
        REJECT="${SAMPLE_OUT}/${base_name}.chrM.cr.rejected.vcf.gz" \
        R="${REF_DIR}/chrM_rCRS.fasta"

    picard MergeVcfs \
        I="${SAMPLE_OUT}/${base_name}.chrM.ncr.reform.vcf.gz" \
        I="${SAMPLE_OUT}/${base_name}.chrM.cr.liftover.vcf.gz" \
        D="${WHOLE_GENOME_REFERENCE_DICT}" \
        O="${SAMPLE_OUT}/${base_name}.chrM.raw.vcf.gz"

    # ---- 最终过滤 ----
    gatk FilterMutectCalls \
        -R "${REF_DIR}/chrM_rCRS.fasta" \
        -V "${SAMPLE_OUT}/${base_name}.chrM.raw.vcf.gz" \
        --stats "${SAMPLE_OUT}/${base_name}.chrM.cr.vcf.gz.stats" \
        --max-alt-allele-count 4 \
        --mitochondria-mode \
        --min-allele-fraction 0.01 \
        -O "${SAMPLE_OUT}/${base_name}.chrM.filtered.vcf.gz"
}

export -f process_sample
export TT_DIR REF_DIR OUT_DIR WHOLE_GENOME_REFERENCE WHOLE_GENOME_REFERENCE_DICT JAVA_OPT THREADS

# 使用 GNU parallel 并行处理
# --colsep 指定分隔符（这里使用制表符或空白字符都可），: ::: 不行，因为我们需要从文件读取
parallel --colsep '\t' process_sample {1} {2} :::: "$LIST_FILE"
