#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONFIG_PATH="${PROJECT_ROOT}/conf/Config.yaml"
MAX_SAMPLES=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --config)
            CONFIG_PATH="$2"
            shift 2
            ;;
        --max-samples)
            MAX_SAMPLES="$2"
            shift 2
            ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

if [[ ! -f "${CONFIG_PATH}" ]]; then
    echo "[ERROR] Config file not found: ${CONFIG_PATH}" >&2
    exit 1
fi

bash "${PROJECT_ROOT}/script/check_env.sh"
source "${PROJECT_ROOT}/script/load_config.sh" "${CONFIG_PATH}"

require_var() {
    local var_name="$1"
    if [[ -z "${!var_name:-}" ]]; then
        echo "[ERROR] Missing required config value: ${var_name}" >&2
        exit 1
    fi
}

for name in \
    CFG_PATHS_LIST_TSV CFG_PATHS_OUTPUT_DIR CFG_PATHS_PYTHON_DIR CFG_PATHS_REFERENCE_FASTA CFG_PATHS_REFERENCE_CACHE_DIR CFG_PATHS_HAPLOGREP_JAR \
    CFG_RUNTIME_CONDA_ENV CFG_RUNTIME_THREADS CFG_RUNTIME_MIN_MAPQ CFG_RUNTIME_MIN_BASEQ CFG_RUNTIME_STRICT_MIN_DP CFG_RUNTIME_STRICT_MIN_VAF CFG_RUNTIME_JAVA_HEAP_GB \
    CFG_TOOLS_CONDA CFG_TOOLS_PYTHON CFG_TOOLS_BWA CFG_TOOLS_SAMTOOLS CFG_TOOLS_BCFTOOLS CFG_TOOLS_JAVA
do
    require_var "${name}"
done

LIST_TSV="${CFG_PATHS_LIST_TSV}"
OUTPUT_DIR="${CFG_PATHS_OUTPUT_DIR}"
PYTHON_DIR="${CFG_PATHS_PYTHON_DIR}"
REFERENCE_FASTA="${CFG_PATHS_REFERENCE_FASTA}"
REFERENCE_CACHE_DIR="${CFG_PATHS_REFERENCE_CACHE_DIR}"
HAPLOGREP_JAR="${CFG_PATHS_HAPLOGREP_JAR}"
CONDA_BIN="${CFG_TOOLS_CONDA}"
CONDA_ENV="${CFG_RUNTIME_CONDA_ENV}"
PYTHON_BIN="${CFG_TOOLS_PYTHON}"
BWA_BIN="${CFG_TOOLS_BWA}"
SAMTOOLS_BIN="${CFG_TOOLS_SAMTOOLS}"
BCFTOOLS_BIN="${CFG_TOOLS_BCFTOOLS}"
JAVA_BIN="${CFG_TOOLS_JAVA}"
THREADS="${CFG_RUNTIME_THREADS}"
MIN_MAPQ="${CFG_RUNTIME_MIN_MAPQ}"
MIN_BASEQ="${CFG_RUNTIME_MIN_BASEQ}"
STRICT_MIN_DP="${CFG_RUNTIME_STRICT_MIN_DP}"
STRICT_MIN_VAF="${CFG_RUNTIME_STRICT_MIN_VAF}"
JAVA_HEAP_GB="${CFG_RUNTIME_JAVA_HEAP_GB}"

if [[ ! -f "${LIST_TSV}" ]]; then
    echo "[ERROR] list.tsv not found: ${LIST_TSV}" >&2
    exit 1
fi

mkdir -p "${OUTPUT_DIR}" "${REFERENCE_CACHE_DIR}"

LOCAL_REFERENCE="${REFERENCE_CACHE_DIR}/chrM_rCRS.fasta"
cp -f "${REFERENCE_FASTA}" "${LOCAL_REFERENCE}"
if [[ -f "${REFERENCE_FASTA}.fai" ]]; then
    cp -f "${REFERENCE_FASTA}.fai" "${LOCAL_REFERENCE}.fai"
fi

if [[ ! -f "${LOCAL_REFERENCE}.bwt" ]]; then
    "${CONDA_BIN}" run -n "${CONDA_ENV}" "${BWA_BIN}" index "${LOCAL_REFERENCE}"
fi
if [[ ! -f "${LOCAL_REFERENCE}.fai" ]]; then
    "${CONDA_BIN}" run -n "${CONDA_ENV}" "${SAMTOOLS_BIN}" faidx "${LOCAL_REFERENCE}"
fi

processed=0
while IFS=$'\t' read -r sample_id fastq_path; do
    if [[ "${sample_id}" == "ID" ]]; then
        continue
    fi
    if [[ -z "${sample_id}" || -z "${fastq_path}" ]]; then
        continue
    fi
    if [[ ! -f "${fastq_path}" ]]; then
        echo "[ERROR] FASTQ not found for ${sample_id}: ${fastq_path}" >&2
        exit 1
    fi

    sample_root="${OUTPUT_DIR}/${sample_id}"
    step1_dir="${sample_root}/1-bam"
    step2_dir="${sample_root}/2-vcf"
    step3_dir="${sample_root}/3-haplogrep3"
    mkdir -p "${step1_dir}" "${step2_dir}" "${step3_dir}"

    rename_map="${step2_dir}/chr_rename.tsv"
    printf "chrM\tMT\n" > "${rename_map}"

    raw_vcf="${step2_dir}/${sample_id}.raw.MT.vcf"
    strict_vcf="${step2_dir}/${sample_id}.strict.MT.vcf"
    raw_haplogrep_txt="${step3_dir}/${sample_id}.raw.haplogrep.txt"
    strict_haplogrep_txt="${step3_dir}/${sample_id}.strict.haplogrep.txt"

    echo "[INFO] Processing ${sample_id}"

    "${CONDA_BIN}" run -n "${CONDA_ENV}" "${PYTHON_BIN}" "${PYTHON_DIR}/align_mt_reads.py" \
        --sample-id "${sample_id}" \
        --fastq "${fastq_path}" \
        --reference "${LOCAL_REFERENCE}" \
        --output-dir "${step1_dir}" \
        --threads "${THREADS}" \
        --bwa-bin "${BWA_BIN}" \
        --samtools-bin "${SAMTOOLS_BIN}" \
        --min-mapq "${MIN_MAPQ}"

    "${CONDA_BIN}" run -n "${CONDA_ENV}" "${PYTHON_BIN}" "${PYTHON_DIR}/call_mt_variants.py" \
        --sample-id "${sample_id}" \
        --input-bam "${step1_dir}/${sample_id}.primary.q${MIN_MAPQ}.bam" \
        --reference "${LOCAL_REFERENCE}" \
        --rename-map "${rename_map}" \
        --output-dir "${step2_dir}" \
        --threads "${THREADS}" \
        --bcftools-bin "${BCFTOOLS_BIN}" \
        --min-mapq "${MIN_MAPQ}" \
        --min-baseq "${MIN_BASEQ}"

    "${CONDA_BIN}" run -n "${CONDA_ENV}" "${PYTHON_BIN}" "${PYTHON_DIR}/filter_mt_vcf_for_haplogrep.py" \
        --input "${raw_vcf}" \
        --output "${strict_vcf}" \
        --stats "${step2_dir}/${sample_id}.strict.filter_stats.tsv" \
        --min-dp "${STRICT_MIN_DP}" \
        --min-vaf "${STRICT_MIN_VAF}" \
        --pass-only

    "${PROJECT_ROOT}/script/run_haplogrep_classify.sh" \
        --conda-bin "${CONDA_BIN}" \
        --conda-env "${CONDA_ENV}" \
        --java-bin "${JAVA_BIN}" \
        --haplogrep-jar "${HAPLOGREP_JAR}" \
        --input-vcf "${raw_vcf}" \
        --output-txt "${raw_haplogrep_txt}" \
        --heap-gb "${JAVA_HEAP_GB}"

    "${PROJECT_ROOT}/script/run_haplogrep_classify.sh" \
        --conda-bin "${CONDA_BIN}" \
        --conda-env "${CONDA_ENV}" \
        --java-bin "${JAVA_BIN}" \
        --haplogrep-jar "${HAPLOGREP_JAR}" \
        --input-vcf "${strict_vcf}" \
        --output-txt "${strict_haplogrep_txt}" \
        --heap-gb "${JAVA_HEAP_GB}"

    "${CONDA_BIN}" run -n "${CONDA_ENV}" "${PYTHON_BIN}" "${PYTHON_DIR}/summarize_mt_results.py" \
        --sample-id "${sample_id}" \
        --raw-vcf "${raw_vcf}" \
        --strict-vcf "${strict_vcf}" \
        --raw-haplogrep "${raw_haplogrep_txt}" \
        --strict-haplogrep "${strict_haplogrep_txt}" \
        --strict-filter-stats "${step2_dir}/${sample_id}.strict.filter_stats.tsv" \
        --coverage-tsv "${step1_dir}/${sample_id}.primary.q${MIN_MAPQ}.coverage.tsv" \
        --output-md "${sample_root}/${sample_id}.summary.md"

    processed=$((processed + 1))
    if [[ -n "${MAX_SAMPLES}" && "${processed}" -ge "${MAX_SAMPLES}" ]]; then
        break
    fi
done < "${LIST_TSV}"

echo "[OK] Pipeline finished for ${processed} sample(s)."
