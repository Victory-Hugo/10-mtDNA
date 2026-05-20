"""
解析 VCF 文件，生成种群等位基因频率矩阵和基因型矩阵，供后续 AMOVA 计算使用。

输出文件：
  allele_freq_matrix.npz  — 频率矩阵（种群 × 位点）和基因型矩阵
  pop_names.txt           — 种群名称列表（按行顺序）
  sample_group_map.tsv    — 过滤后样本 → 种群映射
  vcf_parse_summary.tsv   — 解析统计摘要
"""

import argparse
import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pysam

log = logging.getLogger(__name__)


# ──────────────────────────────────────────────────────────────────────
# 辅助函数
# ──────────────────────────────────────────────────────────────────────

def _read_group_file(group_file: str, id_col: str, pop_col: str) -> pd.DataFrame:
    """读取分组文件，返回仅含 ID 和种群列的 DataFrame（去重、去空）。"""
    path = Path(group_file)
    sep = "\t" if path.suffix in (".tsv", ".txt") else ","
    df = pd.read_csv(path, sep=sep, dtype=str)
    if id_col not in df.columns:
        raise ValueError(f"分组文件缺少列 '{id_col}'，实际列: {list(df.columns)}")
    if pop_col not in df.columns:
        raise ValueError(f"分组文件缺少列 '{pop_col}'，实际列: {list(df.columns)}")
    df = df[[id_col, pop_col]].dropna().drop_duplicates(subset=id_col)
    df.columns = ["sample_id", "pop_name"]
    return df


def _read_vcf_samples(vcf_path: str) -> list[str]:
    """从 VCF 文件头部提取样本 ID 列表。"""
    with pysam.VariantFile(vcf_path, "r") as vcf:
        return list(vcf.header.samples)


def _parse_gt(vcf_path: str, contig_name: str, sample_order: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """
    逐记录读取 VCF，提取基因型矩阵。

    返回:
        gt_matrix   — (n_snps, n_samples) int8，-1=缺失，0=REF，1=ALT
        positions   — (n_snps,) int32，SNP 位置
    """
    # 建立样本名称 → 列索引的映射
    sample_idx = {s: i for i, s in enumerate(sample_order)}
    n_samples = len(sample_order)

    gt_rows: list[np.ndarray] = []
    pos_list: list[int] = []

    with pysam.VariantFile(vcf_path, "r") as vcf:
        for rec in vcf.fetch():
            # 过滤 contig
            if contig_name and rec.chrom != contig_name:
                continue
            gt_row = np.full(n_samples, -1, dtype=np.int8)
            for sname, sidx in sample_idx.items():
                try:
                    sample_gt = rec.samples[sname]["GT"]
                    allele = sample_gt[0]            # 取第一个等位基因（单倍型）
                    if allele is not None:
                        gt_row[sidx] = int(allele)
                except (KeyError, TypeError):
                    pass                             # 保持 -1（缺失）
            gt_rows.append(gt_row)
            pos_list.append(rec.pos)

    if not gt_rows:
        raise ValueError(f"VCF 中未找到 contig='{contig_name}' 的变异记录")

    return np.array(gt_rows, dtype=np.int8), np.array(pos_list, dtype=np.int32)


def _compute_pop_freq(
    gt_matrix: np.ndarray,       # (n_snps, N)
    pop_sizes: np.ndarray,       # (n_pops,)
) -> np.ndarray:
    """
    计算各种群在各位点的 ALT 等位基因频率。

    返回:
        freq_matrix — (n_pops, n_snps) float32，缺失位点频率用全局频率填充
    """
    n_pops = len(pop_sizes)
    n_snps = gt_matrix.shape[0]
    freq = np.zeros((n_pops, n_snps), dtype=np.float32)
    boundaries = np.concatenate([[0], np.cumsum(pop_sizes)]).astype(int)

    for k in range(n_pops):
        sub = gt_matrix[:, boundaries[k]:boundaries[k + 1]]   # (n_snps, n_k)
        valid = sub >= 0
        alt = (sub == 1).sum(axis=1).astype(np.float32)
        n_valid = valid.sum(axis=1).astype(np.float32)
        mask = n_valid > 0
        freq[k, mask] = alt[mask] / n_valid[mask]

    return freq


# ──────────────────────────────────────────────────────────────────────
# 核心业务函数（双模式：可 import 调用）
# ──────────────────────────────────────────────────────────────────────

def run(
    vcf_path: str,
    group_file: str,
    id_col: str,
    pop_col: str,
    contig_name: str,
    min_pop_size: int,
    output_freq_matrix: str,
    output_pop_names: str,
    output_sample_map: str,
    output_summary: str,
) -> dict:
    """
    解析 VCF 并生成等位基因频率矩阵。

    Args:
        vcf_path          — 输入 VCF 文件路径（支持 .gz）
        group_file        — 样本分组 TSV 文件路径
        id_col            — 分组文件中样本 ID 列名
        pop_col           — 分组文件中种群列名
        contig_name       — 仅保留该 CHROM 的位点（空字符串=全部）
        min_pop_size      — 最小种群样本量阈值
        output_freq_matrix — 输出频率矩阵 .npz 文件路径
        output_pop_names  — 输出种群名称 .txt 文件路径
        output_sample_map — 输出样本-种群映射 .tsv 文件路径
        output_summary    — 输出解析摘要 .tsv 文件路径

    Returns:
        统计字典，含 n_samples_vcf、n_pops_retained、n_snps 等字段
    """
    log.info("步骤1：读取分组文件 %s", group_file)
    group_df = _read_group_file(group_file, id_col, pop_col)

    log.info("步骤2：读取 VCF 样本列表 %s", vcf_path)
    vcf_samples = _read_vcf_samples(vcf_path)
    vcf_sample_set = set(vcf_samples)

    # 交集匹配
    matched = group_df[group_df["sample_id"].isin(vcf_sample_set)].copy()
    n_matched = len(matched)
    log.info("VCF 样本数: %d，分组文件样本数: %d，匹配样本数: %d",
             len(vcf_samples), len(group_df), n_matched)

    # 过滤小种群
    pop_counts = matched.groupby("pop_name").size()
    retained_pops = sorted(pop_counts[pop_counts >= min_pop_size].index.tolist())
    if not retained_pops:
        raise ValueError(
            f"过滤后无满足 min_pop_size={min_pop_size} 的种群，"
            f"最大种群大小: {pop_counts.max() if not pop_counts.empty else 0}"
        )
    log.info("种群数（过滤前/后）: %d / %d（min_pop_size=%d）",
             len(pop_counts), len(retained_pops), min_pop_size)

    matched = matched[matched["pop_name"].isin(retained_pops)].copy()

    # 按种群名排序，确定样本顺序
    matched["pop_order"] = matched["pop_name"].map(
        {p: i for i, p in enumerate(retained_pops)}
    )
    matched = matched.sort_values(["pop_order", "sample_id"]).reset_index(drop=True)

    sample_order = matched["sample_id"].tolist()
    pop_sizes = np.array(
        [int((matched["pop_name"] == p).sum()) for p in retained_pops],
        dtype=np.int32,
    )

    log.info("步骤3：解析 VCF 基因型（contig='%s'）…", contig_name)
    gt_matrix, positions = _parse_gt(vcf_path, contig_name, sample_order)
    n_snps_raw = gt_matrix.shape[0]

    # 去除所有样本均缺失的位点
    valid_mask = (gt_matrix >= 0).any(axis=1)
    gt_matrix = gt_matrix[valid_mask]
    positions = positions[valid_mask]
    log.info("SNP 位点数（原始/去除全缺失后）: %d / %d", n_snps_raw, gt_matrix.shape[0])

    log.info("步骤4：计算等位基因频率矩阵…")
    freq_matrix = _compute_pop_freq(gt_matrix, pop_sizes)

    # 创建输出目录
    for out_path in [output_freq_matrix, output_pop_names, output_sample_map, output_summary]:
        Path(out_path).parent.mkdir(parents=True, exist_ok=True)

    # 保存频率矩阵和基因型矩阵（.npz）
    np.savez_compressed(
        output_freq_matrix,
        freq_matrix=freq_matrix,
        gt_matrix=gt_matrix,
        pop_sizes=pop_sizes,
        snp_positions=positions,
    )
    log.info("已保存频率矩阵: %s", output_freq_matrix)

    # 保存种群名称（逐行）
    with open(output_pop_names, "w", encoding="utf-8") as fh:
        fh.write("\n".join(retained_pops) + "\n")
    log.info("已保存种群名称: %s", output_pop_names)

    # 保存样本-种群映射
    matched[["sample_id", "pop_name"]].to_csv(
        output_sample_map, sep="\t", index=False
    )
    log.info("已保存样本映射: %s", output_sample_map)

    # 保存摘要
    summary_data = {
        "n_samples_vcf": len(vcf_samples),
        "n_samples_in_group_file": len(group_df),
        "n_samples_matched": n_matched,
        "n_samples_retained": int(matched.shape[0]),
        "n_pops_before_filter": len(pop_counts),
        "n_pops_retained": len(retained_pops),
        "min_pop_size_threshold": min_pop_size,
        "n_snps_raw": n_snps_raw,
        "n_snps_retained": int(gt_matrix.shape[0]),
        "contig_name": contig_name,
    }
    pd.DataFrame([summary_data]).to_csv(output_summary, sep="\t", index=False)
    log.info("已保存摘要: %s", output_summary)

    log.info("步骤1 完成：%d 种群 × %d SNP", len(retained_pops), gt_matrix.shape[0])
    return summary_data


# ──────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="解析 VCF，生成 AMOVA 所需的等位基因频率矩阵和基因型矩阵"
    )
    p.add_argument("--vcf", required=True, help="输入 VCF 文件路径（支持 .gz）")
    p.add_argument("--group-file", required=True, help="样本分组 TSV 文件路径")
    p.add_argument("--id-col", default="ID", help="分组文件中样本 ID 列名（默认: ID）")
    p.add_argument("--pop-col", required=True, help="分组文件中种群列名")
    p.add_argument("--contig-name", default="", help="仅保留该 CHROM 的位点（空=全部）")
    p.add_argument("--min-pop-size", type=int, default=20,
                   help="最小种群样本量阈值，低于此值的种群被排除（默认: 20）")
    p.add_argument("--output-freq-matrix", required=True,
                   help="输出频率矩阵路径（.npz）")
    p.add_argument("--output-pop-names", required=True,
                   help="输出种群名称列表路径（.txt）")
    p.add_argument("--output-sample-map", required=True,
                   help="输出样本-种群映射路径（.tsv）")
    p.add_argument("--output-summary", required=True,
                   help="输出解析摘要路径（.tsv）")
    return p


def main(argv: list[str] | None = None) -> int:
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )
    args = build_parser().parse_args(argv)
    run(
        vcf_path=args.vcf,
        group_file=args.group_file,
        id_col=args.id_col,
        pop_col=args.pop_col,
        contig_name=args.contig_name,
        min_pop_size=args.min_pop_size,
        output_freq_matrix=args.output_freq_matrix,
        output_pop_names=args.output_pop_names,
        output_sample_map=args.output_sample_map,
        output_summary=args.output_summary,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
