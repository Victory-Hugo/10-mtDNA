#!/usr/bin/env bash
set -uo pipefail

BASEDIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/7-单倍群分型/4-phyloMT"
PHYLOMT="${BASEDIR}/phyloMT.py"
VCF_LIST="${BASEDIR}/conf/vcf.txt"
VCFDIR="${BASEDIR}/output/vcf"
OUTDIR="${BASEDIR}/output"
MERGED_TSV="${OUTDIR}/merged.tsv"
MERGED_EXT="${OUTDIR}/merged.extended.tsv"
# done 状态用单个文本文件记录（每行一个样本名），避免 Windows FS 上大量小文件操作
DONE_FILE="${OUTDIR}/done_samples.txt"
LOGDIR="${VCFDIR}/logs"
TMPDIR="${BASEDIR}/tmp/1-batch-vcf"

TREE="phylotree-new-rcrs@17.3"
HITS=1
JOBS="${JOBS:-16}"
CHUNK_SIZE="${CHUNK_SIZE:-200}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

mkdir -p "${VCFDIR}" "${LOGDIR}" "${TMPDIR}"

command -v parallel >/dev/null 2>&1 || { echo "ERROR: GNU parallel not found" >&2; exit 1; }

if [[ ! -s "${VCF_LIST}" ]]; then
    echo "ERROR: VCF list not found or empty: ${VCF_LIST}" >&2
    exit 1
fi

cd "${BASEDIR}"

# ── Step 1: one-time merge of existing individual TSV files ───────────────────
merge_existing() {
    local header_written=0 count=0
    echo "Merging existing individual TSV files → ${MERGED_TSV} ..."
    while IFS= read -r f; do
        [[ -s "$f" ]] || continue
        if (( header_written == 0 )); then cat "$f"; header_written=1
        else tail -n +2 "$f"
        fi
        (( count++ ))
    done < <(find "${VCFDIR}" -maxdepth 1 -name "*.tsv" ! -name "*.extended.tsv" | sort) \
        >> "${MERGED_TSV}"
    echo "  merged ${count} files."
}

merge_existing_ext() {
    local header_written=0 count=0
    echo "Merging existing extended TSV files → ${MERGED_EXT} ..."
    # 同时收集所有已完成样本名 → DONE_FILE（一次性，比逐文件 touch 快得多）
    > "${DONE_FILE}"
    while IFS= read -r f; do
        [[ -s "$f" ]] || continue
        if (( header_written == 0 )); then cat "$f"; header_written=1
        else tail -n +2 "$f"
        fi
        basename "${f%.extended.tsv}" >> "${DONE_FILE}"
        (( count++ ))
    done < <(find "${VCFDIR}" -maxdepth 1 -name "*.extended.tsv" | sort) \
        >> "${MERGED_EXT}"
    sort -o "${DONE_FILE}" "${DONE_FILE}"
    echo "  merged ${count} files, done list: ${DONE_FILE}"
}

if [[ ! -f "${MERGED_TSV}" ]]; then
    merge_existing
else
    echo "merged.tsv already exists, skipping."
fi

if [[ ! -f "${MERGED_EXT}" ]]; then
    merge_existing_ext
elif [[ ! -f "${DONE_FILE}" ]]; then
    # merged 已存在但 done_samples.txt 缺失：从 extended TSV 文件名重建（一次 find，无 touch）
    echo "Rebuilding done_samples.txt from existing extended TSV filenames..."
    find "${VCFDIR}" -maxdepth 1 -name "*.extended.tsv" -printf "%f\n" \
        | sed 's/\.extended\.tsv$//' \
        | sort > "${DONE_FILE}"
    echo "  $(wc -l < "${DONE_FILE}") entries written."
else
    echo "merged.extended.tsv and done_samples.txt already exist, skipping."
fi

# 确保合并文件有表头（当从空开始时）
if [[ ! -s "${MERGED_TSV}" ]]; then
    printf '"SampleID"\t"Haplogroup"\t"Rank"\t"Quality"\t"Range"\n' > "${MERGED_TSV}"
fi
if [[ ! -s "${MERGED_EXT}" ]]; then
    printf '"SampleID"\t"Haplogroup"\t"Rank"\t"Quality"\t"Range"\t"Input_Sample"\t"Found_Polymorphisms"\t"Missing_Polymorphisms"\t"Extra_Polymorphisms"\t"Found_Weight"\t"Expected_Weight"\t"Sample_Weight"\n' > "${MERGED_EXT}"
fi

# ── Step 2: 构建待处理 VCF 列表（用 comm 差集，避免逐文件检查）────────────────
INPUT_LIST="${TMPDIR}/vcf.all.txt"
awk 'NF && $0 !~ /^[[:space:]]*#/ && $0 ~ /\.vcf(\.gz)?$/ { print $0 }' \
    "${VCF_LIST}" | sort -u > "${INPUT_LIST}"

if [[ ! -s "${INPUT_LIST}" ]]; then
    echo "ERROR: no .vcf / .vcf.gz entries in ${VCF_LIST}" >&2
    exit 1
fi

# 从 VCF 路径提取样本名，与 DONE_FILE 做差集，得到待处理列表
ALL_NAMES="${TMPDIR}/vcf.all_names.txt"
sed 's|.*/||; s/\.vcf\.gz$//; s/\.vcf$//' "${INPUT_LIST}" | sort > "${ALL_NAMES}"

TODO_NAMES="${TMPDIR}/vcf.todo_names.txt"
if [[ -s "${DONE_FILE}" ]]; then
    comm -23 "${ALL_NAMES}" "${DONE_FILE}" > "${TODO_NAMES}"
else
    cp "${ALL_NAMES}" "${TODO_NAMES}"
fi

# 根据 TODO 样本名还原完整 VCF 路径，并过滤不存在的文件
TODO_LIST="${TMPDIR}/vcf.todo.txt"
> "${TODO_LIST}"
# 构建 name→path 映射（awk 一次扫描）
awk 'NR==FNR {
    p=$0; n=p; sub(/.*\//,"",n); sub(/\.vcf\.gz$/,"",n); sub(/\.vcf$/,"",n)
    path[n]=p; next
}
{ if ($0 in path) print path[$0] }' "${INPUT_LIST}" "${TODO_NAMES}" > "${TODO_LIST}"

# 过滤不存在的 VCF
VALID_TODO="${TMPDIR}/vcf.valid_todo.txt"
> "${VALID_TODO}"
while IFS= read -r vcf; do
    if [[ -f "${vcf}" ]]; then
        echo "${vcf}" >> "${VALID_TODO}"
    else
        echo "WARN: VCF not found, skipping: ${vcf}" >&2
    fi
done < "${TODO_LIST}"

TODO_COUNT=$(wc -l < "${VALID_TODO}")
if [[ "${TODO_COUNT}" -eq 0 ]]; then
    echo "All samples already done. Nothing to do."
    exit 0
fi
echo "Remaining: ${TODO_COUNT} VCF files to process."

# ── Step 3: 拆分成 chunk，每个 chunk 一次性加载树 ────────────────────────────
CHUNK_DIR="${TMPDIR}/chunks"
mkdir -p "${CHUNK_DIR}"
rm -f "${CHUNK_DIR}"/chunk_*.txt "${CHUNK_DIR}"/chunk_*.tsv "${CHUNK_DIR}"/chunk_*.extended.tsv

split -l "${CHUNK_SIZE}" -a 4 -d "${VALID_TODO}" "${CHUNK_DIR}/chunk_"
for f in "${CHUNK_DIR}"/chunk_[0-9]*; do
    [[ -f "$f" && "$f" != *.txt ]] && mv "$f" "${f}.txt"
done

CHUNK_LIST="${TMPDIR}/chunks.txt"
find "${CHUNK_DIR}" -name "chunk_*.txt" | sort > "${CHUNK_LIST}"
CHUNK_COUNT=$(wc -l < "${CHUNK_LIST}")
echo "Chunks: ${CHUNK_COUNT} (chunk_size=${CHUNK_SIZE}, jobs=${JOBS})"

# ── Step 4: 并行执行，每个 chunk 调用 phyloMT --input-list ───────────────────
run_one_chunk() {
    local chunk_list="$1"
    local chunk_name
    chunk_name="$(basename "${chunk_list%.txt}")"
    local out_tsv="${CHUNK_DIR}/${chunk_name}.tsv"

    "${PYTHON_BIN}" "${PHYLOMT}" \
        --tree "${TREE}" \
        --input-list "${chunk_list}" \
        --output "${out_tsv}" \
        --hits "${HITS}" \
        --extended-report \
        > "${LOGDIR}/${chunk_name}.stdout.log" \
        2> "${LOGDIR}/${chunk_name}.stderr.log"
}

export BASEDIR PHYLOMT CHUNK_DIR LOGDIR TREE HITS PYTHON_BIN
export -f run_one_chunk

parallel \
    --jobs "${JOBS}" \
    --halt never \
    --joblog "${LOGDIR}/parallel.$(date +%Y%m%d_%H%M%S).joblog" \
    run_one_chunk :::: "${CHUNK_LIST}"

# ── Step 5: 合并 chunk 结果，更新 done 记录 ───────────────────────────────────
echo "Merging chunk results..."
while IFS= read -r chunk_list; do
    chunk_name="$(basename "${chunk_list%.txt}")"
    out_tsv="${CHUNK_DIR}/${chunk_name}.tsv"
    out_ext="${CHUNK_DIR}/${chunk_name}.extended.tsv"

    if [[ ! -s "${out_tsv}" ]]; then
        echo "WARN: chunk ${chunk_name} produced no output, skipping." >&2
        continue
    fi

    tail -n +2 "${out_tsv}" >> "${MERGED_TSV}"
    [[ -s "${out_ext}" ]] && tail -n +2 "${out_ext}" >> "${MERGED_EXT}"

    # 追加完成的样本名到 DONE_FILE（一次写，不用逐文件 touch）
    sed 's|.*/||; s/\.vcf\.gz$//; s/\.vcf$//' "${chunk_list}" >> "${DONE_FILE}"

    rm -f "${out_tsv}" "${out_ext}" "${chunk_list}"
done < "${CHUNK_LIST}"

# 对 DONE_FILE 重新排序去重，保持查找效率
sort -u -o "${DONE_FILE}" "${DONE_FILE}"

echo "Done. Merged outputs: ${MERGED_TSV}  ${MERGED_EXT}"
