#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONFIG_PATH="${PROJECT_ROOT}/conf/Config.yaml"

source "${SCRIPT_DIR}/load_config.sh" "${CONFIG_PATH}"

FASTQ="${1:-${CFG_PATHS_FASTQ_DIR}/SRR36389141.fastq.gz}"
SAMPLE_OUTDIR="${2:-${CFG_PATHS_OUTPUT_DIR}/SRR36389141}"
THREADS="${THREADS:-${CFG_RUNTIME_THREADS}}"

REF_SOURCE="${CFG_PATHS_REFERENCE_FASTA}"
HAPLOGREP_JAR="${CFG_PATHS_HAPLOGREP_JAR}"
HAPLOGREP_DIR="$(dirname "${HAPLOGREP_JAR}")"
PY_FILTER="${CFG_PATHS_PYTHON_DIR}/filter_mt_vcf_for_haplogrep.py"

mkdir -p "${SAMPLE_OUTDIR}"

sample_name="$(basename "${FASTQ}")"
sample_name="${sample_name%.fastq.gz}"
sample_name="${sample_name%.fq.gz}"
sample_name="${sample_name%.fastq}"
sample_name="${sample_name%.fq}"

LOCAL_REF="${SAMPLE_OUTDIR}/chrM_rCRS.fasta"
LOCAL_REF_FAI="${SAMPLE_OUTDIR}/chrM_rCRS.fasta.fai"

cp -f "${REF_SOURCE}" "${LOCAL_REF}"
if [[ -f "${REF_SOURCE}.fai" ]]; then
    cp -f "${REF_SOURCE}.fai" "${LOCAL_REF_FAI}"
else
    "${CFG_TOOLS_CONDA}" run -n "${CFG_RUNTIME_CONDA_ENV}" "${CFG_TOOLS_SAMTOOLS}" faidx "${LOCAL_REF}"
fi

if [[ ! -f "${LOCAL_REF}.bwt" ]]; then
    "${CFG_TOOLS_CONDA}" run -n "${CFG_RUNTIME_CONDA_ENV}" "${CFG_TOOLS_BWA}" index "${LOCAL_REF}"
fi

rename_map="${SAMPLE_OUTDIR}/chr_rename.tsv"
printf "chrM\tMT\n" > "${rename_map}"

aligned_bam="${SAMPLE_OUTDIR}/${sample_name}.rCRS.sorted.bam"
primary_bam="${SAMPLE_OUTDIR}/${sample_name}.rCRS.primary.q20.bam"
raw_vcfgz="${SAMPLE_OUTDIR}/${sample_name}.bcftools.raw.chrM.vcf.gz"
raw_mt_vcf="${SAMPLE_OUTDIR}/${sample_name}.bcftools.raw.MT.vcf"
strict_vcf="${SAMPLE_OUTDIR}/${sample_name}.bcftools.strict.MT.vcf"
strict_stats="${SAMPLE_OUTDIR}/${sample_name}.bcftools.strict.filter_stats.tsv"
haplo_raw_txt="${SAMPLE_OUTDIR}/${sample_name}.bcftools.raw.haplogrep.txt"
haplo_strict_txt="${SAMPLE_OUTDIR}/${sample_name}.bcftools.strict.haplogrep.txt"

"${CFG_TOOLS_CONDA}" run -n "${CFG_RUNTIME_CONDA_ENV}" bash -lc "
set -euo pipefail
${CFG_TOOLS_BWA} mem -t ${THREADS} -M \
  -R '@RG\tID:${sample_name}\tSM:${sample_name}\tPL:ILLUMINA\tLB:lib1\tPU:unit1' \
  '${LOCAL_REF}' '${FASTQ}' \
| ${CFG_TOOLS_SAMTOOLS} sort -@ ${THREADS} -o '${aligned_bam}' -
${CFG_TOOLS_SAMTOOLS} index '${aligned_bam}'
${CFG_TOOLS_SAMTOOLS} flagstat '${aligned_bam}' > '${SAMPLE_OUTDIR}/${sample_name}.rCRS.flagstat.txt'
${CFG_TOOLS_SAMTOOLS} idxstats '${aligned_bam}' > '${SAMPLE_OUTDIR}/${sample_name}.rCRS.idxstats.txt'
${CFG_TOOLS_SAMTOOLS} view -@ ${THREADS} -b -F 2308 -q ${CFG_RUNTIME_MIN_MAPQ} '${aligned_bam}' -o '${primary_bam}'
${CFG_TOOLS_SAMTOOLS} index '${primary_bam}'
${CFG_TOOLS_SAMTOOLS} flagstat '${primary_bam}' > '${SAMPLE_OUTDIR}/${sample_name}.rCRS.primary.q20.flagstat.txt'
${CFG_TOOLS_SAMTOOLS} coverage '${primary_bam}' > '${SAMPLE_OUTDIR}/${sample_name}.rCRS.primary.q20.coverage.tsv'
${CFG_TOOLS_BCFTOOLS} mpileup -f '${LOCAL_REF}' -a FORMAT/DP,FORMAT/AD -q ${CFG_RUNTIME_MIN_MAPQ} -Q ${CFG_RUNTIME_MIN_BASEQ} -Ou '${primary_bam}' \
| ${CFG_TOOLS_BCFTOOLS} call -mv --ploidy 1 -Ou \
| ${CFG_TOOLS_BCFTOOLS} norm -f '${LOCAL_REF}' -m -both -Oz -o '${raw_vcfgz}'
${CFG_TOOLS_BCFTOOLS} index -f '${raw_vcfgz}'
${CFG_TOOLS_BCFTOOLS} annotate --rename-chrs '${rename_map}' '${raw_vcfgz}' -Ov -o '${raw_mt_vcf}'
${CFG_TOOLS_BCFTOOLS} stats '${raw_vcfgz}' > '${SAMPLE_OUTDIR}/${sample_name}.bcftools.raw.stats.txt'
"

"${CFG_TOOLS_CONDA}" run -n "${CFG_RUNTIME_CONDA_ENV}" "${CFG_TOOLS_PYTHON}" "${PY_FILTER}" \
    --input "${raw_mt_vcf}" \
    --output "${strict_vcf}" \
    --stats "${strict_stats}" \
    --min-dp "${CFG_RUNTIME_STRICT_MIN_DP}" \
    --min-vaf "${CFG_RUNTIME_STRICT_MIN_VAF}" \
    --pass-only

"${CFG_TOOLS_CONDA}" run -n "${CFG_RUNTIME_CONDA_ENV}" bash -lc "
set -euo pipefail
cd '${HAPLOGREP_DIR}'
${CFG_TOOLS_JAVA} -Xmx${CFG_RUNTIME_JAVA_HEAP_GB}G -jar '${HAPLOGREP_JAR}' classify \
  --input '${raw_mt_vcf}' \
  --tree phylotree-rcrs@17.2 \
  --output '${haplo_raw_txt}' \
  --extend-report \
  --write-fasta
${CFG_TOOLS_JAVA} -Xmx${CFG_RUNTIME_JAVA_HEAP_GB}G -jar '${HAPLOGREP_JAR}' classify \
  --input '${strict_vcf}' \
  --tree phylotree-rcrs@17.2 \
  --output '${haplo_strict_txt}' \
  --extend-report \
  --write-fasta
"
