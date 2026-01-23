#!/usr/bin/env bash
set -euo pipefail

input_dir="/mnt/f/onedrive/文档（科研）/脚本/Download/10-mtDNA/5-经典群体遗传学计算/output/group_Population"
output_pdf="${input_dir}/genetics_summary.pdf"
script_dir="/mnt/f/onedrive/文档（科研）/脚本/Download/10-mtDNA/5-经典群体遗传学计算/python/vis.R"

Rscript "${script_dir}" \
    "${input_dir}" \
    "${output_pdf}"
