#!/bin/bash
#? phyloMT 使用示例
#? 注意！fasta文件可以使用任何坐标系（rCRS/RSRS），但 HSD/VCF 文件必须与树的坐标系一致，否则会导致分型错误。
# 树选择原则：
#   phylotree-new-rcrs@17.2  — 新版 JSON 树，rCRS 坐标系（适用于 rCRS 坐标的 FASTA/VCF/HSD）
#   phylotree-new-rsrs@17.2  — 新版 JSON 树，RSRS 坐标系（适用于 RSRS 坐标的 FASTA/VCF/HSD）
#   phylotree-rcrs@17.2      — 旧版 XML 树，rCRS 坐标系
#   phylotree-rsrs@17.1      — 旧版 XML 树，RSRS 坐标系

BASEDIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/7-单倍群分型/4-phyloMT"
PHYLOMT="${BASEDIR}/phyloMT.py"
THREADS=16
HITS=1
cd "${BASEDIR}"

# ── 示例 1：FASTA 输入 × 新版 rCRS 树 ──────────────────────────────────────
python3 "${PHYLOMT}" \
    --tree phylotree-new-rcrs@17.3 \
    --input "${BASEDIR}/input/fasta/example.fasta" \
    --output "${BASEDIR}/output/新版/17.3/fasta/example_fasta_rcrs_output.tsv" \
    --threads ${THREADS} --hits ${HITS} --extended-report

# ── 示例 2：FASTA 输入 × 新版 RSRS 树 ──────────────────────────────────────
python3 "${PHYLOMT}" \
    --tree phylotree-new-rsrs@17.3 \
    --input "${BASEDIR}/input/fasta/example.fasta" \
    --output "${BASEDIR}/output/新版/17.3/fasta/example_fasta_rsrs_output.tsv" \
    --threads ${THREADS} --hits ${HITS} --extended-report

# ── 示例 3：HSD 输入 × 新版 rCRS 树（HSD 变异为 rCRS 坐标）────────────────
python3 "${PHYLOMT}" \
    --tree phylotree-new-rcrs@17.3 \
    --input "${BASEDIR}/input/hsd/example.hsd" \
    --output "${BASEDIR}/output/新版/17.3/hsd/example_hsd_rcrs_output.tsv" \
    --threads ${THREADS} --hits ${HITS} --extended-report



# ── 示例 4：WGS VCF 输入 × 新版 rCRS 树（VCF 以 rCRS 为参考）──────────────
python3 "${PHYLOMT}" \
    --tree phylotree-new-rcrs@17.3 \
    --input "${BASEDIR}/input/vcf/example-wgs.vcf" \
    --output "${BASEDIR}/output/新版/17.3/vcf/example_wgs_rcrs_output.tsv" \
    --threads ${THREADS} --hits ${HITS} --extended-report



# ── 示例 5：芯片 VCF 输入 × 新版 rCRS 树（--chip 模式）─────────────────────
python3 "${PHYLOMT}" \
    --tree phylotree-new-rcrs@17.3 \
    --input "${BASEDIR}/input/vcf/example-microarray.vcf" \
    --output "${BASEDIR}/output/新版/17.3/vcf/example_microarray_rcrs_output.tsv" \
    --threads ${THREADS} --hits ${HITS} --extended-report --chip


# ── 示例 6：FASTA 输入 × 旧版 rCRS 树 ──────────────────────────────────────
# python3 "${PHYLOMT}" \
#     --tree phylotree-rcrs@17.2 \
#     --input "${BASEDIR}/input/fasta/example.fasta" \
#     --output "${BASEDIR}/output/旧版/fasta/example_fasta_rcrs_output.tsv" \
#     --threads ${THREADS} --hits ${HITS} --extended-report

# ── 示例 7：HSD 输入 × 旧版 rCRS 树 ────────────────────────────────────────
# python3 "${PHYLOMT}" \
#     --tree phylotree-rcrs@17.2 \
#     --input "${BASEDIR}/input/hsd/example.hsd" \
#     --output "${BASEDIR}/output/旧版/hsd/example_hsd_rcrs_output.tsv" \
#     --threads ${THREADS} --hits ${HITS} --extended-report

# ── 示例 8：WGS VCF 输入 × 旧版 rCRS 树 ────────────────────────────────────
# python3 "${PHYLOMT}" \
#     --tree phylotree-rcrs@17.2 \
#     --input "${BASEDIR}/input/vcf/example-wgs.vcf" \
#     --output "${BASEDIR}/output/旧版/vcf/example_wgs_rcrs_output.tsv" \
#     --threads ${THREADS} --hits ${HITS} --extended-report
