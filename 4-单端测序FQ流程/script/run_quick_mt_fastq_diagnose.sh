#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONFIG_PATH="${PROJECT_ROOT}/conf/Config.yaml"

source "${SCRIPT_DIR}/load_config.sh" "${CONFIG_PATH}"

FASTQ="${1:-${CFG_PATHS_FASTQ_DIR}/SRR36389141.fastq.gz}"
OUTDIR="${2:-${CFG_PATHS_OUTPUT_DIR}}"
REFERENCE="${3:-${CFG_PATHS_SHIFTED_REFERENCE_FASTA}}"
THREADS="${THREADS:-${CFG_RUNTIME_THREADS}}"

"${CFG_TOOLS_CONDA}" run -n "${CFG_RUNTIME_CONDA_ENV}" "${CFG_TOOLS_PYTHON}" "${CFG_PATHS_PYTHON_DIR}/quick_mt_fastq_diagnose.py" \
    --fastq "${FASTQ}" \
    --outdir "${OUTDIR}" \
    --reference "${REFERENCE}" \
    --threads "${THREADS}"
