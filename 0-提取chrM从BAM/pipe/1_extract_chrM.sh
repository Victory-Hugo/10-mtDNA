#!/bin/bash
set -euo pipefail

# ==================================================
# 基础路径定义
BASE_DIR="/mnt/f/Onedrive/文档（科研）/脚本/Download/10-mtDNA/4-Mutect2"
LIST_FILE="${BASE_DIR}/meta/list.txt"
OUT_DIR="${BASE_DIR}/output"
LOG_DIR="${OUT_DIR}/logs"
SUCCESS_LOG="${LOG_DIR}/extract_chrM.success.log"
THREADS=8
# ==================================================

mkdir -p "${LOG_DIR}"
touch "${SUCCESS_LOG}"

declare -A SUCCESS_SAMPLES
while IFS=$'\t' read -r sample_id _rest; do
    [[ -z "${sample_id}" ]] && continue
    SUCCESS_SAMPLES["${sample_id}"]=1
done < "${SUCCESS_LOG}"

process_sample() {
    local base_name bam_file
    base_name="$1"
    bam_file="$2"

    if [[ -n "${SUCCESS_SAMPLES["${base_name}"]+x}" ]]; then
        echo "Skipping sample (already success): ${base_name}"
        return 0
    fi

    echo "Processing sample: ${base_name}"
    mkdir -p "${OUT_DIR}/${base_name}"

    samtools view -@ "${THREADS}" -b \
        -o "${OUT_DIR}/${base_name}/${base_name}.chrM.bam" \
        "${bam_file}" chrM

    samtools index -@ "${THREADS}" \
        "${OUT_DIR}/${base_name}/${base_name}.chrM.bam"

    {
        flock -x 9
        if ! grep -Fxq "${base_name}" "${SUCCESS_LOG}"; then
            echo "${base_name}" >> "${SUCCESS_LOG}"
        fi
    } 9>>"${SUCCESS_LOG}"
}

export -f process_sample
export OUT_DIR LOG_DIR SUCCESS_LOG THREADS

while IFS=$'\t' read -r sample_id bam_path _rest; do
    [[ -z "${sample_id}" ]] && continue
    process_sample "${sample_id}" "${bam_path}"
done < "${LIST_FILE}"
