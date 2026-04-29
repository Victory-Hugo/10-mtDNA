#!/usr/bin/env bash
# Benchmark: old (per-sample) vs new (--input-list batch) mode
set -euo pipefail

BASEDIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/7-单倍群分型/4-phyloMT"
PHYLOMT="${BASEDIR}/phyloMT.py"
PYTHON_BIN="${PYTHON_BIN:-python3}"
TREE="phylotree-new-rcrs@17.3"
HITS=1

WGS_VCF="${BASEDIR}/input/vcf/example-wgs.vcf"
SPLIT_DIR="${BASEDIR}/tmp/benchmark_split"
OUTDIR="${BASEDIR}/output/新版/17.3"
OLD_OUTDIR="${OUTDIR}/vcf/benchmark_old"
NEW_OUTDIR="${OUTDIR}/vcf/benchmark_new"

mkdir -p "${SPLIT_DIR}" "${OLD_OUTDIR}" "${NEW_OUTDIR}"

# ── Step 1: 将 WGS VCF 拆分为单样本 VCF ──────────────────────────────────────
echo "=== Step 1: 拆分多样本 VCF（bcftools）==="
SAMPLE_NAMES=$("${PYTHON_BIN}" -c "
import pysam
vf = pysam.VariantFile('${WGS_VCF}')
print('\n'.join(vf.header.samples))
")
SAMPLE_COUNT=$(echo "${SAMPLE_NAMES}" | wc -l)
echo "  样本数: ${SAMPLE_COUNT}"

echo "${SAMPLE_NAMES}" | while IFS= read -r sample; do
    out="${SPLIT_DIR}/${sample}.vcf"
    [[ -f "${out}" ]] && continue
    bcftools view -s "${sample}" -O v -o "${out}" "${WGS_VCF}" 2>/dev/null
done
echo "  拆分完成。"

# 生成拆分文件列表
SPLIT_LIST="${BASEDIR}/tmp/benchmark_split_list.txt"
find "${SPLIT_DIR}" -name "*.vcf" | sort > "${SPLIT_LIST}"
SAMPLE_COUNT=$(wc -l < "${SPLIT_LIST}")
echo "  共 ${SAMPLE_COUNT} 个单样本 VCF"

# ── Step 2: 旧模式 benchmark（逐个调用，顺序执行）────────────────────────────
echo ""
echo "=== Step 2: 旧模式（每次 1 个 VCF，共 ${SAMPLE_COUNT} 次调用）==="
rm -f "${OLD_OUTDIR}"/*.tsv

OLD_START=$(date +%s%3N)
while IFS= read -r vcf; do
    name="$(basename "${vcf%.vcf}")"
    "${PYTHON_BIN}" "${PHYLOMT}" \
        --tree "${TREE}" \
        --input "${vcf}" \
        --output "${OLD_OUTDIR}/${name}.tsv" \
        --hits "${HITS}" \
        --extended-report \
        > /dev/null 2>&1
done < "${SPLIT_LIST}"
OLD_END=$(date +%s%3N)
OLD_MS=$(( OLD_END - OLD_START ))

echo "  旧模式耗时: ${OLD_MS} ms ($(( OLD_MS / 1000 ))s)"
OLD_ROWS=$(cat "${OLD_OUTDIR}"/*.tsv | grep -v '^"SampleID"' | wc -l)
echo "  旧模式结果行数: ${OLD_ROWS}"

# ── Step 3: 新模式 benchmark（--input-list，一次调用）───────────────────────
echo ""
echo "=== Step 3: 新模式（--input-list，一次调用处理 ${SAMPLE_COUNT} 个 VCF）==="
NEW_OUT="${NEW_OUTDIR}/batch_output.tsv"
rm -f "${NEW_OUTDIR}"/*.tsv

NEW_START=$(date +%s%3N)
"${PYTHON_BIN}" "${PHYLOMT}" \
    --tree "${TREE}" \
    --input-list "${SPLIT_LIST}" \
    --output "${NEW_OUT}" \
    --hits "${HITS}" \
    --extended-report \
    > /dev/null 2>&1
NEW_END=$(date +%s%3N)
NEW_MS=$(( NEW_END - NEW_START ))

echo "  新模式耗时: ${NEW_MS} ms ($(( NEW_MS / 1000 ))s)"
NEW_ROWS=$(tail -n +2 "${NEW_OUT}" | wc -l)
echo "  新模式结果行数: ${NEW_ROWS}"

# ── Step 4: 对比输出是否一致 ────────────────────────────────────────────────
echo ""
echo "=== Step 4: 对比输出一致性 ==="
"${PYTHON_BIN}" - <<'PYEOF'
import os, glob

old_dir = "/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/7-单倍群分型/4-phyloMT/output/新版/17.3/vcf/benchmark_old"
new_file = "/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/7-单倍群分型/4-phyloMT/output/新版/17.3/vcf/benchmark_new/batch_output.tsv"
ref_file = "/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/7-单倍群分型/4-phyloMT/output/新版/17.3/vcf/example_wgs_rcrs_output.tsv"

def load_tsv(path):
    rows = {}
    with open(path) as f:
        lines = f.readlines()
    for line in lines[1:]:
        line = line.strip().strip('"')
        if not line:
            continue
        cols = [c.strip('"') for c in line.split('\t')]
        if cols:
            rows[cols[0]] = cols  # key = SampleID
    return rows

# 从旧模式所有独立文件合并
old_rows = {}
for tsv in sorted(glob.glob(os.path.join(old_dir, "*.tsv"))):
    if ".extended." in tsv:
        continue
    old_rows.update(load_tsv(tsv))

new_rows = load_tsv(new_file)
ref_rows = load_tsv(ref_file) if os.path.exists(ref_file) else {}

print(f"  旧模式样本数: {len(old_rows)}")
print(f"  新模式样本数: {len(new_rows)}")
if ref_rows:
    print(f"  参考输出样本数: {len(ref_rows)}")

# 逐样本对比旧 vs 新
mismatches = []
for sid in sorted(old_rows):
    if sid not in new_rows:
        mismatches.append(f"  MISSING in new: {sid}")
        continue
    o = old_rows[sid]
    n = new_rows[sid]
    # 比较 Haplogroup(1) 和 Quality(3)
    if o[1] != n[1] or o[3] != n[3]:
        mismatches.append(f"  MISMATCH {sid}: old={o[1]}/{o[3]} new={n[1]}/{n[3]}")

if mismatches:
    print(f"  !! {len(mismatches)} 个差异:")
    for m in mismatches[:10]:
        print(m)
else:
    print("  旧模式 vs 新模式: 完全一致 ✓")

# 新模式 vs 参考输出
if ref_rows:
    ref_mismatches = []
    for sid in sorted(ref_rows):
        if sid not in new_rows:
            ref_mismatches.append(f"  MISSING in new: {sid}")
            continue
        r = ref_rows[sid]
        n = new_rows[sid]
        if r[1] != n[1] or r[3] != n[3]:
            ref_mismatches.append(f"  MISMATCH {sid}: ref={r[1]}/{r[3]} new={n[1]}/{n[3]}")
    if ref_mismatches:
        print(f"  !! 新模式 vs 参考: {len(ref_mismatches)} 个差异:")
        for m in ref_mismatches[:10]:
            print(m)
    else:
        print("  新模式 vs 参考输出: 完全一致 ✓")
PYEOF

# ── Step 5: 汇总 ─────────────────────────────────────────────────────────────
echo ""
echo "=== 性能汇总 ==="
printf "  旧模式: %6d ms (%d 次树加载)\n" "${OLD_MS}" "${SAMPLE_COUNT}"
printf "  新模式: %6d ms (%d 次树加载)\n" "${NEW_MS}" 1
if (( NEW_MS > 0 )); then
    SPEEDUP=$(python3 -c "print(f'{${OLD_MS}/${NEW_MS}:.1f}x')")
    echo "  加速比: ${SPEEDUP}"
fi
