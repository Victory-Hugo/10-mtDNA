#!/bin/bash

set -euo pipefail


INPUT_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/7-单倍群分型/1-HaploCart/example/1.fa"
sample_prefix="1"
OUTPUT_FILE="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/7-单倍群分型/1-HaploCart/output/HaploCart"
HC_FILE_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/7-单倍群分型/1-HaploCart/data"
posterior_file="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/7-单倍群分型/1-HaploCart/output/HaploCart_posterior"
THREADS="16"
VGAN="vgan"

"$VGAN" haplocart \
    --hc-files "$HC_FILE_DIR" \
    -f "$INPUT_FILE" \
    -o "$OUTPUT_FILE" \
    -pf "$posterior_file" \
    -t "$THREADS" \
    -s "$sample_prefix"
    