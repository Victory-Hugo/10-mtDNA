#!/usr/bin/env bash
set -euo pipefail

JAVA_BIN=""
CONDA_BIN=""
CONDA_ENV=""
HAPLOGREP_JAR=""
INPUT_VCF=""
OUTPUT_TXT=""
TREE_NAME="phylotree-rcrs@17.2"
HEAP_GB="4"

while [[ $# -gt 0 ]]; do
    case "$1" in
        --java-bin)
            JAVA_BIN="$2"
            shift 2
            ;;
        --conda-bin)
            CONDA_BIN="$2"
            shift 2
            ;;
        --conda-env)
            CONDA_ENV="$2"
            shift 2
            ;;
        --haplogrep-jar)
            HAPLOGREP_JAR="$2"
            shift 2
            ;;
        --input-vcf)
            INPUT_VCF="$2"
            shift 2
            ;;
        --output-txt)
            OUTPUT_TXT="$2"
            shift 2
            ;;
        --tree-name)
            TREE_NAME="$2"
            shift 2
            ;;
        --heap-gb)
            HEAP_GB="$2"
            shift 2
            ;;
        *)
            echo "[ERROR] Unknown argument: $1" >&2
            exit 1
            ;;
    esac
done

for required_var in JAVA_BIN HAPLOGREP_JAR INPUT_VCF OUTPUT_TXT TREE_NAME HEAP_GB; do
    if [[ -z "${!required_var}" ]]; then
        echo "[ERROR] Missing required argument: ${required_var}" >&2
        exit 1
    fi
done

HAPLOGREP_DIR="$(cd "$(dirname "${HAPLOGREP_JAR}")" && pwd)"
mkdir -p "$(dirname "${OUTPUT_TXT}")"

(
    cd "${HAPLOGREP_DIR}"
    if [[ -n "${CONDA_ENV}" ]]; then
        "${CONDA_BIN}" run -n "${CONDA_ENV}" "${JAVA_BIN}" -Xmx"${HEAP_GB}"G -jar "${HAPLOGREP_JAR}" classify \
            --input "${INPUT_VCF}" \
            --tree "${TREE_NAME}" \
            --output "${OUTPUT_TXT}" \
            --extend-report \
            --write-fasta
    else
        "${JAVA_BIN}" -Xmx"${HEAP_GB}"G -jar "${HAPLOGREP_JAR}" classify \
            --input "${INPUT_VCF}" \
            --tree "${TREE_NAME}" \
            --output "${OUTPUT_TXT}" \
            --extend-report \
            --write-fasta
    fi
)
