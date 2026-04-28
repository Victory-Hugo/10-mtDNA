#!/usr/bin/env bash
set -uo pipefail

BASEDIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/7-单倍群分型/4-phyloMT"
PHYLOMT="${BASEDIR}/phyloMT.py"
VCF_LIST="${BASEDIR}/conf/vcf.txt"
VCFDIR="${BASEDIR}/output/vcf"       # legacy individual TSVs + done markers
OUTDIR="${BASEDIR}/output"
MERGED_TSV="${OUTDIR}/merged.tsv"
MERGED_EXT="${OUTDIR}/merged.extended.tsv"
DONEDIR="${VCFDIR}/.done"
LOGDIR="${VCFDIR}/logs"
TMPDIR="${BASEDIR}/tmp/1-batch-vcf"

TREE="phylotree-new-rcrs@17.3"
THREADS=1
HITS=1
JOBS="${JOBS:-16}"
PYTHON_BIN="${PYTHON_BIN:-python3}"

mkdir -p "${VCFDIR}" "${DONEDIR}" "${LOGDIR}" "${TMPDIR}"

command -v parallel >/dev/null 2>&1 || { echo "ERROR: GNU parallel not found" >&2; exit 1; }
command -v flock    >/dev/null 2>&1 || { echo "ERROR: flock not found" >&2; exit 1; }

if [[ ! -s "${VCF_LIST}" ]]; then
    echo "ERROR: VCF list not found or empty: ${VCF_LIST}" >&2
    exit 1
fi

cd "${BASEDIR}"

# ── Step 1: one-time merge of existing individual TSV files ───────────────────
merge_existing() {
    local header_written=0
    local count=0
    echo "Merging existing individual TSV files → ${MERGED_TSV} ..."
    while IFS= read -r f; do
        [[ -s "$f" ]] || continue
        if (( header_written == 0 )); then
            cat "$f" >> "${MERGED_TSV}"
            header_written=1
        else
            tail -n +2 "$f" >> "${MERGED_TSV}"
        fi
        (( count++ ))
    done < <(find "${VCFDIR}" -maxdepth 1 -name "*.tsv" ! -name "*.extended.tsv" | sort)
    echo "  merged ${count} files."
}

merge_existing_ext() {
    local header_written=0
    local count=0
    echo "Merging existing extended TSV files → ${MERGED_EXT} ..."
    while IFS= read -r f; do
        [[ -s "$f" ]] || continue
        if (( header_written == 0 )); then
            cat "$f" >> "${MERGED_EXT}"
            header_written=1
        else
            tail -n +2 "$f" >> "${MERGED_EXT}"
        fi
        # create done marker so these samples are skipped later
        local name
        name="$(basename "${f%.extended.tsv}")"
        touch "${DONEDIR}/${name}.done"
        (( count++ ))
    done < <(find "${VCFDIR}" -maxdepth 1 -name "*.extended.tsv" | sort)
    echo "  merged ${count} files."
}

if [[ ! -f "${MERGED_TSV}" ]]; then
    merge_existing
else
    echo "merged.tsv already exists, skipping merge step."
fi

if [[ ! -f "${MERGED_EXT}" ]]; then
    merge_existing_ext
else
    echo "merged.extended.tsv already exists, skipping merge step."
    # still populate done markers if missing
    find "${VCFDIR}" -maxdepth 1 -name "*.extended.tsv" | while IFS= read -r f; do
        local_name="$(basename "${f%.extended.tsv}")"
        touch "${DONEDIR}/${local_name}.done"
    done
fi

# ── Step 2: build header lines (needed when appending to empty merged files) ──
# Write header to merged files if they were just created empty
if [[ ! -s "${MERGED_TSV}" ]]; then
    printf '"SampleID"\t"Haplogroup"\t"Rank"\t"Quality"\t"Range"\n' > "${MERGED_TSV}"
fi
if [[ ! -s "${MERGED_EXT}" ]]; then
    printf '"SampleID"\t"Haplogroup"\t"Rank"\t"Quality"\t"Range"\t"Input_Sample"\t"Found_Polymorphisms"\t"Missing_Polymorphisms"\t"Extra_Polymorphisms"\t"Found_Weight"\t"Expected_Weight"\t"Sample_Weight"\n' > "${MERGED_EXT}"
fi

# ── Step 3: prepare input list ────────────────────────────────────────────────
INPUT_LIST="${TMPDIR}/vcf.inputs.txt"
awk 'NF && $0 !~ /^[[:space:]]*#/ && $0 ~ /\.vcf(\.gz)?$/ { print $0 }' \
    "${VCF_LIST}" | sort -u > "${INPUT_LIST}"

if [[ ! -s "${INPUT_LIST}" ]]; then
    echo "ERROR: no .vcf / .vcf.gz entries in ${VCF_LIST}" >&2
    exit 1
fi

# ── Step 4: per-sample function ───────────────────────────────────────────────
run_one_vcf() {
    local vcf="$1"
    local name tmp_tsv tmp_ext

    name="$(basename "${vcf}")"
    name="${name%.vcf.gz}"
    name="${name%.vcf}"

    # skip if already done
    if [[ -f "${DONEDIR}/${name}.done" ]]; then
        echo "SKIP ${name}"
        return 0
    fi

    # skip if VCF missing (warn, don't abort)
    if [[ ! -f "${vcf}" ]]; then
        echo "WARN: VCF not found, skipping: ${vcf}" >&2
        return 0
    fi

    tmp_tsv="${TMPDIR}/${name}.tsv"
    tmp_ext="${TMPDIR}/${name}.extended.tsv"

    # run phyloMT; on failure log and continue
    if ! "${PYTHON_BIN}" "${PHYLOMT}" \
            --tree "${TREE}" \
            --input "${vcf}" \
            --output "${tmp_tsv}" \
            --threads "${THREADS}" \
            --hits "${HITS}" \
            --extended-report \
            > "${LOGDIR}/${name}.stdout.log" \
            2> "${LOGDIR}/${name}.stderr.log"; then
        echo "ERROR: phyloMT failed for ${name} (see ${LOGDIR}/${name}.stderr.log)" >&2
        return 0
    fi

    # append data rows (skip header) to merged files under a lock
    if [[ -s "${tmp_tsv}" ]]; then
        (flock 200; tail -n +2 "${tmp_tsv}" >> "${MERGED_TSV}") 200>"${MERGED_TSV}.lock"
    fi
    if [[ -s "${tmp_ext}" ]]; then
        (flock 201; tail -n +2 "${tmp_ext}" >> "${MERGED_EXT}") 201>"${MERGED_EXT}.lock"
    fi

    # mark done and clean up tmp
    touch "${DONEDIR}/${name}.done"
    rm -f "${tmp_tsv}" "${tmp_ext}"
}

export BASEDIR PHYLOMT DONEDIR LOGDIR TMPDIR MERGED_TSV MERGED_EXT TREE THREADS HITS PYTHON_BIN
export -f run_one_vcf

# ── Step 5: run in parallel, continue on errors ───────────────────────────────
parallel \
    --jobs "${JOBS}" \
    --halt never \
    --joblog "${LOGDIR}/parallel.$(date +%Y%m%d_%H%M%S).joblog" \
    run_one_vcf :::: "${INPUT_LIST}"

echo "Done. Merged outputs: ${MERGED_TSV}  ${MERGED_EXT}"
