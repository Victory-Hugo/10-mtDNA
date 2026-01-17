#!/bin/bash
set -euo pipefail

# ==================================================
# 基础路径与参数定义（从配置文件读取）
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONF_FILE="${SCRIPT_DIR}/../conf/2_Mutect2.conf"
if [[ ! -f "${CONF_FILE}" ]]; then
    echo "Config not found: ${CONF_FILE}" >&2
    exit 1
fi
source "${CONF_FILE}"
# ==================================================

# # 如果需要建立参考基因组索引（只需执行一次），可以取消下面注释
# samtools faidx "${REF_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta"
# ${PICARD} CreateSequenceDictionary \
#     R="${WHOLE_GENOME_REFERENCE}" \
#     O="${WHOLE_GENOME_REFERENCE_DICT}"

mkdir -p "${LOG_DIR}"
touch "${SUCCESS_LOG}"

declare -A SUCCESS_SAMPLES
while IFS=$'\t' read -r sample_id _rest; do
    [[ -z "${sample_id}" ]] && continue
    SUCCESS_SAMPLES["${sample_id}"]=1
done < "${SUCCESS_LOG}"

# 单个样本处理函数
process_sample() {
    local base_name bam_file cram_index
    base_name="$1"
    bam_file="$2"

    if [[ -n "${SUCCESS_SAMPLES["${base_name}"]+x}" ]]; then
        echo "Skipping sample (already success): ${base_name}"
        return 0
    fi

    echo "Processing sample: ${base_name}"

    # 为当前样本创建输出目录
    SAMPLE_OUT="${OUT_DIR}/${base_name}"
    mkdir -p "${SAMPLE_OUT}"

    # ---- 先抽取chrM再补RG（避免对全基因组BAM做RG）----
    if [[ "${bam_file}" == *.cram ]]; then
        cram_index="${bam_file}.crai"
        if [[ ! -f "${cram_index}" ]]; then
            samtools index -@ ${THREADS} "${bam_file}"
        fi
    fi
    mkdir -p "${TMP_DIR}"
    TMP_CHRM="${TMP_DIR}/${base_name}.chrM.bam"
    samtools view -@ ${THREADS} -b -T "${WHOLE_GENOME_REFERENCE}" \
        "${bam_file}" chrM -o "${TMP_CHRM}"
    TMP_BAM="${TMP_DIR}/${base_name}.chrM.rg.bam"
    samtools addreplacerg -w -@ ${THREADS} \
        -r "ID:${base_name}\tSM:${base_name}\tPL:illumina\tLB:lib1\tPU:unit1" \
        -o "${TMP_BAM}" "${TMP_CHRM}"

    # ---- 原始数据处理 ----
    # TMP_BAM 已只含 chrM，避免 -L 触发索引要求
    gatk --java-options "${JAVA_OPT}" PrintReads \
        -I "${TMP_BAM}" \
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter \
        -O "${SAMPLE_OUT}/${base_name}.chrM.bam"

    samtools sort -@ ${THREADS} -O bam \
        -o "${SAMPLE_OUT}/${base_name}.chrM.sorted.bam" \
        "${SAMPLE_OUT}/${base_name}.chrM.bam"

    ${PICARD} BuildBamIndex \
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

    # 补齐contig长度，避免与liftover后的VCF字典不一致
    gatk UpdateVCFSequenceDictionary \
        -V "${SAMPLE_OUT}/${base_name}.chrM.ncr.reform.vcf.gz" \
        -O "${SAMPLE_OUT}/${base_name}.chrM.ncr.reform.dict.vcf.gz" \
        --sequence-dictionary "${REF_DIR}/chrM_rCRS.dict" \
        --replace true

    # ---- 数据重处理流程 ----
    gatk RevertSam \
        -I "${SAMPLE_OUT}/${base_name}.chrM.sorted.bam" \
        -O "${SAMPLE_OUT}/${base_name}.chrM.reverted.bam"

    ${PICARD} SamToFastq \
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
    ${PICARD} AddOrReplaceReadGroups \
        I="${SAMPLE_OUT}/${base_name}.chrM.shifted.sorted.bam" \
        O="${SAMPLE_OUT}/${base_name}.chrM.shifted.rg.bam" \
        RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM="${base_name}"

    ${PICARD} MarkDuplicates \
        INPUT="${SAMPLE_OUT}/${base_name}.chrM.shifted.rg.bam" \
        OUTPUT="${SAMPLE_OUT}/${base_name}.chrM.dedup.bam" \
        METRICS_FILE="${SAMPLE_OUT}/${base_name}.chrM.dedup.metrics"

    ${PICARD} BuildBamIndex \
        INPUT="${SAMPLE_OUT}/${base_name}.chrM.dedup.bam" \
        OUTPUT="${SAMPLE_OUT}/${base_name}.chrM.dedup.bai"

    ${PICARD} CollectWgsMetrics \
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
    ${PICARD} LiftoverVcf \
        I="${SAMPLE_OUT}/${base_name}.chrM.cr.vcf.gz" \
        O="${SAMPLE_OUT}/${base_name}.chrM.cr.liftover.vcf.gz" \
        CHAIN="${REF_DIR}/chrM_rCRS.shift_back.chain" \
        REJECT="${SAMPLE_OUT}/${base_name}.chrM.cr.rejected.vcf.gz" \
        R="${REF_DIR}/chrM_rCRS.fasta"

    ${PICARD} MergeVcfs \
        I="${SAMPLE_OUT}/${base_name}.chrM.ncr.reform.dict.vcf.gz" \
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

    {
        flock -x 9
        if ! grep -Fxq "${base_name}" "${SUCCESS_LOG}"; then
            echo "${base_name}" >> "${SUCCESS_LOG}"
        fi
    } 9>>"${SUCCESS_LOG}"

    # ---- 成功后清理临时BAM ----
    rm -f "${TMP_BAM}" "${TMP_BAM}.bai" "${TMP_CHRM}"
}

export -f process_sample
export REF_DIR OUT_DIR WHOLE_GENOME_REFERENCE WHOLE_GENOME_REFERENCE_DICT JAVA_OPT THREADS LOG_DIR SUCCESS_LOG PICARD TMP_DIR


# 使用 GNU parallel 并行处理
# --colsep 指定分隔符（这里使用制表符或空白字符都可），: ::: 不行，因为我们需要从文件读取
parallel --jobs "${PARALLEL_JOBS}" --colsep '\t' process_sample {1} {2} :::: "$LIST_FILE"
