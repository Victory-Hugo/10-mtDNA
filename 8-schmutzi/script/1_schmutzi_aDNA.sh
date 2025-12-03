#!/usr/bin/bash
# batch_schmutzi.sh
# Usage: ./batch_schmutzi.sh sample_list.tsv OUTPUT_DIR TEMP_DIR BASE_DIR
# sample_list.tsv 需包含两列：<BAM路径><TAB><SampleID>；输出和临时目录需手动指定；BASE_DIR 由调用脚本传入

set -euo pipefail

############################################
# 全局配置（从调用脚本的环境变量读取）
############################################
THREADS=${THREADS:-8}
LENGTH_DEAM=${LENGTH_DEAM:-5}
LIB_TYPE=${LIB_TYPE:-"double"}
MIN_READS=${MIN_READS:-50}
USE_ENDO=${USE_ENDO:-1}
ITERATIONS=${ITERATIONS:-2}
LOG2FASTA_QUAL=${LOG2FASTA_QUAL:-20}
LOG2FASTA_INDEL=${LOG2FASTA_INDEL:-0}
RESUME=${RESUME:-0}
############################################

if [[ $# -ne 4 ]]; then
  echo "Usage: $0 <sample_list.tsv> <output_dir> <temp_dir> <base_dir>" >&2
  exit 1
fi

LIST_FILE="$1"
OUTPUT_ROOT="$2"
TEMP_ROOT="$3"
ROOT="$4"

if [[ ! -s $LIST_FILE ]]; then
  echo "Error: input list '$LIST_FILE' not found or empty." >&2
  exit 1
fi

mkdir -p "$OUTPUT_ROOT" "$TEMP_ROOT"
OUTPUT_ROOT="$(realpath "$OUTPUT_ROOT")"
TEMP_ROOT="$(realpath "$TEMP_ROOT")"
ROOT="$(realpath "$ROOT")"
REF=${REF:-"$ROOT/share/schmutzi/refs/human_MT.fa"}
AL_FREQ_DIR=${AL_FREQ_DIR:-"$ROOT/share/schmutzi/alleleFreqMT/eurasian/freqs"}
SCHMUTZI_PERL=${SCHMUTZI_PERL:-"$ROOT/src/schmutzi.pl"}
LOG2FASTA=${LOG2FASTA:-"$ROOT/src/log2fasta"}
SAMTOOLS_PATH=${SAMTOOLS_PATH:-"samtools"}
BWA_PATH=${BWA_PATH:-"bwa"}

# 准备参考序列长度
if [[ ! -f "${REF}.fai" ]]; then
  "$SAMTOOLS_PATH" faidx "$REF"
fi
REF_NAME=$(cut -f1 "$REF.fai" | head -1)
REF_LENGTH=$(awk 'NR==1{print $2}' "$REF.fai")

#--------------------------------------------------------------------
process_sample() {
  local bam_path="$1"
  local sample_id="$2"

  local SAMPLE_OUTPUT_DIR="${OUTPUT_ROOT}/${sample_id}"
  local FINAL_FASTA="${SAMPLE_OUTPUT_DIR}/${sample_id}.fasta"

  # 断点续跑检查：如果输出文件已存在且非空，则跳过
  if [[ "$RESUME" -eq 1 && -s "$FINAL_FASTA" ]]; then
    echo "[$(date +%F\ %T)] SKIP: ${sample_id} (已完成, 输出文件已存在)"
    return
  fi

  bam_path="$(realpath "$bam_path")"
  if [[ ! -f $bam_path ]]; then
    echo "Error: BAM file '$bam_path' not found." >&2
    exit 1
  fi

  local INPUT_DIR="$ROOT/data"
  local SAMPLE_TEMP_DIR="${TEMP_ROOT}/${sample_id}"
  local SAMPLE_OUTPUT_DIR="${OUTPUT_ROOT}/${sample_id}"

  mkdir -p "$INPUT_DIR" "$SAMPLE_TEMP_DIR" "$SAMPLE_OUTPUT_DIR"

  # 复制 BAM
  local DEST_BAM="$INPUT_DIR/${sample_id}.bam"
  if [[ "$(realpath "$bam_path")" != "$(realpath "$DEST_BAM")" ]]; then
    cp "$bam_path" "$DEST_BAM"
  fi

  # 中间文件路径
  local FASTQ_FILE="$SAMPLE_TEMP_DIR/${sample_id}.fastq"
  local SAM_FILE="$SAMPLE_TEMP_DIR/${sample_id}.sam"
  local RAW_BAM="$SAMPLE_TEMP_DIR/${sample_id}.raw.bam"
  local SORTED_BAM="$SAMPLE_TEMP_DIR/${sample_id}.sorted.bam"
  local CALMD_BAM="$SAMPLE_TEMP_DIR/${sample_id}.calmd.bam"

  echo "[$(date +%F\ %T)] Processing sample: ${sample_id}"

  # 1. BAM → FASTQ
  "$SAMTOOLS_PATH" fastq "$DEST_BAM" > "$FASTQ_FILE"

  # 2. 比对
  "$BWA_PATH" mem -t "$THREADS" "$REF" "$FASTQ_FILE" > "$SAM_FILE"

  # 3. SAM → BAM → sort → calmd
  "$SAMTOOLS_PATH" view -bS "$SAM_FILE" > "$RAW_BAM"
  "$SAMTOOLS_PATH" sort "$RAW_BAM" -o "$SORTED_BAM"
  "$SAMTOOLS_PATH" index "$SORTED_BAM"
  "$SAMTOOLS_PATH" calmd -b "$SORTED_BAM" "$REF" > "$CALMD_BAM"
  "$SAMTOOLS_PATH" index "$CALMD_BAM"

  # 4. 检查覆盖度
  local num_reads
  num_reads=$("$SAMTOOLS_PATH" view -c "$CALMD_BAM")
  if [[ "$num_reads" -lt "$MIN_READS" ]]; then
    echo "[WARNING] 样本 ${sample_id} 只有 ${num_reads} 条 reads (<${MIN_READS})，输出全N序列"
    # 生成长度与参考相同的 N 序列 FASTA
    local lowcov_fa="$SAMPLE_OUTPUT_DIR/${sample_id}.fasta"
    {
      echo ">${sample_id}"
      printf 'N%.0s' $(seq 1 "$REF_LENGTH")
      echo
    } > "$lowcov_fa"
    echo "[$(date +%F\ %T)] Completed low-coverage: ${sample_id}"
    return
  fi

  # 5. contDeam：捕获失败
  local NOENDO_OPT=""
  if [[ "$USE_ENDO" -eq 0 ]]; then
    NOENDO_OPT="--noendo"
  fi
  echo "[DEBUG] Running contDeam with: lengthDeam=$LENGTH_DEAM, library=$LIB_TYPE, use_endo=$USE_ENDO"
  if ! "$ROOT/src/contDeam.pl" \
      --lengthDeam "$LENGTH_DEAM" --library "$LIB_TYPE" $NOENDO_OPT \
      --out "${SAMPLE_TEMP_DIR}/${sample_id}_deam" "$REF" "$CALMD_BAM" 2>&1 | tee "${SAMPLE_TEMP_DIR}/${sample_id}_contDeam.log"; then
    echo "[WARNING] contDeam 失败，详见 ${SAMPLE_TEMP_DIR}/${sample_id}_contDeam.log"
    local failfa="$SAMPLE_OUTPUT_DIR/${sample_id}.fasta"
    {
      echo ">${sample_id}"
      printf 'N%.0s' $(seq 1 "$REF_LENGTH")
      echo
    } > "$failfa"
    echo "[$(date +%F\ %T)] Completed contDeam-fail: ${sample_id}"
    return
  fi

  # 6. schmutzi
  "$SCHMUTZI_PERL" \
    --iterations "$ITERATIONS" \
    --t "$THREADS" \
    --uselength \
    --ref "$REF" \
    --out "${SAMPLE_TEMP_DIR}/${sample_id}_wpred" \
    "${SAMPLE_TEMP_DIR}/${sample_id}_deam" \
    "$AL_FREQ_DIR" \
    "$CALMD_BAM"
    # 如需关闭污染预测，添加 --notusepredC

  # 7. log2fasta & 结果整理
  "$LOG2FASTA" \
    -name "$sample_id" \
    -q "$LOG2FASTA_QUAL" -indel "$LOG2FASTA_INDEL" \
    "${SAMPLE_TEMP_DIR}/${sample_id}_wpred_final_endo.log" \
    > "${SAMPLE_TEMP_DIR}/${sample_id}_wpred_final_endo.HQ.fa"

  mv "${SAMPLE_TEMP_DIR}/${sample_id}_wpred_final_endo.HQ.fa" \
     "${SAMPLE_OUTPUT_DIR}/${sample_id}.fasta"
  mv "${SAMPLE_TEMP_DIR}/${sample_id}_wpred_final_endo.log" \
     "${SAMPLE_OUTPUT_DIR}/"

  echo "[$(date +%F\ %T)] Completed: ${sample_id}"
}
#--------------------------------------------------------------------

while IFS=$'\t' read -r bam_path sample_id; do
  [[ -z "${bam_path}" || "${bam_path:0:1}" == "#" ]] && continue
  if [[ -z "${sample_id}" ]]; then
    echo "Error: line with BAM '$bam_path' lacks a second column (sample ID)." >&2
    exit 1
  fi
  process_sample "$bam_path" "$sample_id"
done < "$LIST_FILE"
