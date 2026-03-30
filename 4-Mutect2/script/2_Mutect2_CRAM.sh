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

preflight_check() {
    local errors=0

    # 检查软件可用性
    for tool in "${SAMTOOLS}" "${GATK}" "${BWA}" "${BGZIP}" "${BCFTOOLS}" "${PARALLEL}"; do
        if ! command -v "${tool}" &>/dev/null; then
            echo "[preflight] ERROR: 工具不可用: ${tool}" >&2
            errors=$(( errors + 1 ))
        fi
    done

    # 检查 picard jar
    local picard_jar
    picard_jar=$(echo "${PICARD}" | awk '{for(i=1;i<=NF;i++){if($i=="-jar"){print $(i+1)}}}')
    if [[ -n "${picard_jar}" && ! -f "${picard_jar}" ]]; then
        echo "[preflight] ERROR: picard.jar 不存在: ${picard_jar}" >&2
        errors=$(( errors + 1 ))
    fi

    # 检查参考序列文件
    for ref_file in \
        "${REF_CHRM}" "${REF_CHRM}.fai" "${REF_CHRM_DICT}" \
        "${REF_CHRM_SHIFTED}" "${REF_CHRM_SHIFTED}.fai" "${REF_CHRM_SHIFTED_DICT}" \
        "${REF_SHIFT_CHAIN}" \
        "${WHOLE_GENOME_REFERENCE}" "${WHOLE_GENOME_REFERENCE_DICT}"; do
        if [[ ! -f "${ref_file}" ]]; then
            echo "[preflight] ERROR: 参考文件不存在: ${ref_file}" >&2
            errors=$(( errors + 1 ))
        fi
    done

    # 检查样本列表文件
    if [[ ! -f "${LIST_FILE}" ]]; then
        echo "[preflight] ERROR: 样本列表不存在: ${LIST_FILE}" >&2
        errors=$(( errors + 1 ))
    else
        local sample_count
        sample_count=$(wc -l < "${LIST_FILE}")
        echo "[preflight] 样本列表: ${LIST_FILE} (${sample_count} 行)"
    fi

    if (( errors > 0 )); then
        echo "[preflight] 发现 ${errors} 个错误，退出。" >&2
        exit 1
    fi
    echo "[preflight] 所有检查通过。"
}

# ==================================================

# # 如果需要建立参考基因组索引（只需执行一次），可以取消下面注释
# ${SAMTOOLS} faidx "${WHOLE_GENOME_REFERENCE}"
# ${PICARD} CreateSequenceDictionary \
#     R="${WHOLE_GENOME_REFERENCE}" \
#     O="${WHOLE_GENOME_REFERENCE_DICT}"

preflight_check

mkdir -p "${LOG_DIR}"
touch "${SUCCESS_LOG}"
touch "${FAIL_LOG}"

declare -A SUCCESS_SAMPLES
while IFS=$'\t' read -r field1 field2; do
    [[ -z "${field1}" ]] && continue
    # 兼容新旧格式：旧格式只有样本名（field2为空），新格式为"时间戳\t样本名"
    if [[ -n "${field2}" ]]; then
        SUCCESS_SAMPLES["${field2}"]=1
    else
        SUCCESS_SAMPLES["${field1}"]=1
    fi
done < "${SUCCESS_LOG}"

# 单个样本处理函数
process_sample() {
    set -euo pipefail
    local base_name bam_file cram_index
    local _step=""
    base_name="$1"
    bam_file="$2"

    trap 'echo -e "$(date +"%Y-%m-%d %H:%M:%S")\t${base_name}\t${_step:-unknown}" >> "${FAIL_LOG}"' ERR

    if [[ -n "${SUCCESS_SAMPLES["${base_name}"]+x}" ]]; then
        echo "Skipping sample (already success): ${base_name}"
        return 0
    fi

    echo "Processing sample: ${base_name}"

    # 为当前样本创建输出目录
    SAMPLE_OUT="${OUT_DIR}/${base_name}"
    mkdir -p "${SAMPLE_OUT}"

    # ---- 先抽取chrM再补RG（避免对全基因组BAM做RG）----
    _step="index_cram"
    if [[ "${bam_file}" == *.cram ]]; then
        cram_index="${bam_file}.crai"
        if [[ ! -f "${cram_index}" ]]; then
            ${SAMTOOLS} index -@ ${THREADS} "${bam_file}"
        fi
    fi
    mkdir -p "${TMP_DIR}"
    TMP_CHRM="${TMP_DIR}/${base_name}.chrM.bam"
    _step="extract_chrM"
    ${SAMTOOLS} view -@ ${THREADS} -b -T "${WHOLE_GENOME_REFERENCE}" \
        "${bam_file}" chrM -o "${TMP_CHRM}"
    TMP_BAM="${TMP_DIR}/${base_name}.chrM.rg.bam"
    _step="add_rg"
    ${SAMTOOLS} addreplacerg -w -@ ${THREADS} \
        -r "ID:${base_name}\tSM:${base_name}\tPL:illumina\tLB:lib1\tPU:unit1" \
        -o "${TMP_BAM}" "${TMP_CHRM}"

    # ---- 原始数据处理 ----
    # TMP_BAM 已只含 chrM，避免 -L 触发索引要求
    _step="PrintReads"
    ${GATK} --java-options "${JAVA_OPT}" PrintReads \
        -I "${TMP_BAM}" \
        --read-filter MateOnSameContigOrNoMappedMateReadFilter \
        --read-filter MateUnmappedAndUnmappedReadFilter \
        -O "${SAMPLE_OUT}/${base_name}.chrM.bam"

    _step="sort_bam"
    ${SAMTOOLS} sort -@ ${THREADS} -O bam \
        -o "${SAMPLE_OUT}/${base_name}.chrM.sorted.bam" \
        "${SAMPLE_OUT}/${base_name}.chrM.bam"

    _step="BuildBamIndex_sorted"
    ${PICARD} BuildBamIndex \
        INPUT="${SAMPLE_OUT}/${base_name}.chrM.sorted.bam" \
        OUTPUT="${SAMPLE_OUT}/${base_name}.chrM.sorted.bai"

    # ---- 非控制区变异检测 ----
    _step="Mutect2_ncr"
    ${GATK} --java-options "${JAVA_OPT}" Mutect2 \
        -R "${REF_CHRM}" \
        -L chrM \
        --mitochondria-mode \
        --annotation StrandBiasBySample \
        --max-reads-per-alignment-start 75 \
        --max-mnp-distance 0 \
        -I "${SAMPLE_OUT}/${base_name}.chrM.sorted.bam" \
        -O "${SAMPLE_OUT}/${base_name}.chrM.ncr.vcf"

    _step="sed_ncr"
    sed '/^##contig.*\(random\|Un\|_alt\|HLA\|EBV\|chr[1-9XY]\)/d' \
        "${SAMPLE_OUT}/${base_name}.chrM.ncr.vcf" \
        > "${SAMPLE_OUT}/${base_name}.chrM.ncr.reform.vcf"

    _step="bgzip_ncr"
    ${BGZIP} -f -@ ${THREADS} "${SAMPLE_OUT}/${base_name}.chrM.ncr.reform.vcf"
    _step="bcftools_index_ncr"
    ${BCFTOOLS} index --threads ${THREADS} -f -t \
        "${SAMPLE_OUT}/${base_name}.chrM.ncr.reform.vcf.gz"

    # 补齐contig长度，避免与liftover后的VCF字典不一致
    _step="UpdateVCFSequenceDictionary"
    ${GATK} UpdateVCFSequenceDictionary \
        -V "${SAMPLE_OUT}/${base_name}.chrM.ncr.reform.vcf.gz" \
        -O "${SAMPLE_OUT}/${base_name}.chrM.ncr.reform.dict.vcf.gz" \
        --sequence-dictionary "${REF_CHRM_DICT}" \
        --replace true

    # ---- 数据重处理流程 ----
    _step="RevertSam"
    ${GATK} RevertSam \
        -I "${SAMPLE_OUT}/${base_name}.chrM.sorted.bam" \
        -O "${SAMPLE_OUT}/${base_name}.chrM.reverted.bam"

    _step="SamToFastq"
    ${PICARD} SamToFastq \
        I="${SAMPLE_OUT}/${base_name}.chrM.reverted.bam" \
        FASTQ="${SAMPLE_OUT}/${base_name}_R1.chrM.fq.gz" \
        SECOND_END_FASTQ="${SAMPLE_OUT}/${base_name}_R2.chrM.fq.gz"

    # ---- 重比对流程 ----
    _step="bwa_mem"
    ${BWA} mem -t ${THREADS} -M \
        -R "@RG\tID:group1\tSM:${base_name}\tPL:illumina\tLB:CasCADE\tPU:unit1" \
        "${REF_CHRM_SHIFTED}" \
        "${SAMPLE_OUT}/${base_name}_R1.chrM.fq.gz" \
        "${SAMPLE_OUT}/${base_name}_R2.chrM.fq.gz" \
        > "${SAMPLE_OUT}/${base_name}.chrM.shifted.sam"

    _step="sort_shifted"
    ${SAMTOOLS} sort -@ ${THREADS} -O bam \
        -o "${SAMPLE_OUT}/${base_name}.chrM.shifted.sorted.bam" \
        "${SAMPLE_OUT}/${base_name}.chrM.shifted.sam"

    # ---- 数据标记和指标收集 ----
    _step="AddOrReplaceReadGroups"
    ${PICARD} AddOrReplaceReadGroups \
        I="${SAMPLE_OUT}/${base_name}.chrM.shifted.sorted.bam" \
        O="${SAMPLE_OUT}/${base_name}.chrM.shifted.rg.bam" \
        RGID=4 RGLB=lib1 RGPL=illumina RGPU=unit1 RGSM="${base_name}"

    _step="MarkDuplicates"
    ${PICARD} MarkDuplicates \
        INPUT="${SAMPLE_OUT}/${base_name}.chrM.shifted.rg.bam" \
        OUTPUT="${SAMPLE_OUT}/${base_name}.chrM.dedup.bam" \
        METRICS_FILE="${SAMPLE_OUT}/${base_name}.chrM.dedup.metrics"

    _step="BuildBamIndex_dedup"
    ${PICARD} BuildBamIndex \
        INPUT="${SAMPLE_OUT}/${base_name}.chrM.dedup.bam" \
        OUTPUT="${SAMPLE_OUT}/${base_name}.chrM.dedup.bai"

    _step="CollectWgsMetrics"
    ${PICARD} CollectWgsMetrics \
        I="${SAMPLE_OUT}/${base_name}.chrM.dedup.bam" \
        O="${SAMPLE_OUT}/${base_name}.chrM.wgs_metrics.txt" \
        R="${REF_CHRM_SHIFTED}"

    # ---- 控制区变异检测 ----
    _step="Mutect2_cr"
    ${GATK} Mutect2 \
        -R "${REF_CHRM_SHIFTED}" \
        -L chrM \
        --mitochondria-mode \
        --annotation StrandBiasBySample \
        --max-reads-per-alignment-start 75 \
        --max-mnp-distance 0 \
        -I "${SAMPLE_OUT}/${base_name}.chrM.dedup.bam" \
        -O "${SAMPLE_OUT}/${base_name}.chrM.cr.vcf.gz"

    # ---- 变异数据整合 ----
    _step="LiftoverVcf"
    ${PICARD} LiftoverVcf \
        I="${SAMPLE_OUT}/${base_name}.chrM.cr.vcf.gz" \
        O="${SAMPLE_OUT}/${base_name}.chrM.cr.liftover.vcf.gz" \
        CHAIN="${REF_SHIFT_CHAIN}" \
        REJECT="${SAMPLE_OUT}/${base_name}.chrM.cr.rejected.vcf.gz" \
        R="${REF_CHRM}"

    _step="MergeVcfs"
    ${PICARD} MergeVcfs \
        I="${SAMPLE_OUT}/${base_name}.chrM.ncr.reform.dict.vcf.gz" \
        I="${SAMPLE_OUT}/${base_name}.chrM.cr.liftover.vcf.gz" \
        D="${WHOLE_GENOME_REFERENCE_DICT}" \
        O="${SAMPLE_OUT}/${base_name}.chrM.raw.vcf.gz"

    # ---- 最终过滤 ----
    _step="FilterMutectCalls"
    ${GATK} FilterMutectCalls \
        -R "${REF_CHRM}" \
        -V "${SAMPLE_OUT}/${base_name}.chrM.raw.vcf.gz" \
        --stats "${SAMPLE_OUT}/${base_name}.chrM.cr.vcf.gz.stats" \
        --max-alt-allele-count 4 \
        --mitochondria-mode \
        --min-allele-fraction 0.01 \
        -O "${SAMPLE_OUT}/${base_name}.chrM.filtered.vcf.gz"

    {
        flock -x 9
        if ! awk -F'\t' -v s="${base_name}" '{if($NF==s){found=1}} END{exit !found}' "${SUCCESS_LOG}" 2>/dev/null; then
            printf '%s\t%s\n' "$(date '+%Y-%m-%d %H:%M:%S')" "${base_name}" >> "${SUCCESS_LOG}"
        fi
    } 9>>"${SUCCESS_LOG}"

    # ---- 成功后清理临时BAM ----
    rm -f "${TMP_BAM}" "${TMP_BAM}.bai" "${TMP_CHRM}"
}

export -f process_sample
export REF_DIR OUT_DIR WHOLE_GENOME_REFERENCE WHOLE_GENOME_REFERENCE_DICT JAVA_OPT THREADS LOG_DIR SUCCESS_LOG PICARD TMP_DIR
export SAMTOOLS GATK BWA BGZIP BCFTOOLS FAIL_LOG REF_CHRM REF_CHRM_DICT REF_CHRM_SHIFTED REF_CHRM_SHIFTED_DICT REF_SHIFT_CHAIN


# 使用 GNU parallel 并行处理
# --colsep 指定分隔符（这里使用制表符或空白字符都可），: ::: 不行，因为我们需要从文件读取
parallel --jobs "${PARALLEL_JOBS}" --colsep '\t' process_sample {1} {2} :::: "$LIST_FILE"
