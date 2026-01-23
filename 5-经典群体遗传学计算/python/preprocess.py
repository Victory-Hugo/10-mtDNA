from __future__ import annotations

from pathlib import Path
from typing import Dict, Iterable, List, Tuple
import logging

import numpy as np
import pandas as pd
import allel

logger = logging.getLogger(__name__)


def read_sample_table(table_path: Path) -> pd.DataFrame:
    logger.info("读取样本信息表: %s", table_path)
    if table_path.suffix.lower() == ".tsv":
        df = pd.read_csv(table_path, sep="\t")
        logger.info("样本表读取完成: %s 行, %s 列", df.shape[0], df.shape[1])
        return df
    if table_path.suffix.lower() == ".csv":
        df = pd.read_csv(table_path)
        logger.info("样本表读取完成: %s 行, %s 列", df.shape[0], df.shape[1])
        return df
    raise ValueError("样本信息表格仅支持 .csv 或 .tsv 格式")


def validate_sample_table(sample_df: pd.DataFrame, id_col: str, group_cols: Iterable[str]) -> None:
    logger.info("校验样本表字段: ID 列=%s, 分群列=%s", id_col, list(group_cols))
    missing_cols = [col for col in [id_col, *group_cols] if col not in sample_df.columns]
    if missing_cols:
        raise ValueError(f"样本信息表格缺少列: {missing_cols}")


def summarize_groups(sample_df: pd.DataFrame, group_col: str) -> pd.DataFrame:
    logger.info("汇总分群列: %s", group_col)
    group_counts = (
        sample_df[[group_col]]
        .dropna()
        .value_counts()
        .reset_index()
        .rename(columns={group_col: "group", 0: "count"})
    )
    logger.info("分群列 %s 可用群体数: %s", group_col, len(group_counts))
    return group_counts


def check_group_counts(group_counts: pd.DataFrame, group_col: str) -> None:
    if group_counts.empty:
        raise ValueError(f"分群列 {group_col} 没有可用的群体信息")

    min_count = group_counts["count"].min()
    all_ge_20 = (group_counts["count"] >= 20).all()

    logger.info("分群列 %s 的群体数量: %s", group_col, len(group_counts))
    logger.info("分群列 %s 统计结果:\n%s", group_col, group_counts.to_string(index=False))

    if min_count < 2:
        raise ValueError(f"分群列 {group_col} 中存在样本数 < 2 的群体，已终止")

    if all_ge_20:
        logger.info("分群列 %s: 所有群体样本数 >= 20，统计充分，继续计算", group_col)
    else:
        logger.warning("分群列 %s: 所有群体样本数 >= 2，将继续计算（样本量可能偏小）", group_col)


def read_vcf(vcf_path: Path) -> Dict[str, np.ndarray]:
    logger.info("读取 VCF: %s", vcf_path)
    if vcf_path.suffix.lower() not in {".vcf", ".gz"}:
        raise ValueError("VCF 文件需为 .vcf 或 .vcf.gz")

    try:
        vcf = allel.read_vcf(
            str(vcf_path),
            fields=[
                "samples",
                "variants/CHROM",
                "variants/POS",
                "variants/REF",
                "variants/ALT",
                "variants/numalt",
                "calldata/GT",
            ],
        )
        logger.info("VCF 读取完成: 样本数 %s, 位点数 %s", len(vcf["samples"]), len(vcf["variants/POS"]))
        return vcf
    except Exception as exc:
        raise ValueError(f"VCF 读取失败: {exc}") from exc


def analyze_variant_types(vcf: Dict[str, np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
    logger.debug("分析变异类型 (InDel/多等位)")
    ref = vcf["variants/REF"].astype(str)
    alt = vcf["variants/ALT"]
    numalt = vcf["variants/numalt"]

    ref_len = np.char.str_len(ref)
    alt_len = np.vectorize(lambda x: len(x) if x is not None else 0)(alt)
    has_indel = (ref_len != 1) | (alt_len > 1).any(axis=1)
    is_multiallelic = numalt > 1

    return has_indel, is_multiallelic


def filter_variants(
    vcf: Dict[str, np.ndarray],
    chrom_name: str,
    skip_complex_variants: bool,
) -> Dict[str, np.ndarray]:
    chrom = vcf["variants/CHROM"].astype(str)
    if chrom_name not in np.unique(chrom):
        raise ValueError(f"VCF 不包含染色体 {chrom_name}")

    chrom_mask = chrom == chrom_name

    has_indel, is_multiallelic = analyze_variant_types(vcf)

    if has_indel.any():
        logger.warning("VCF 包含 InDel 变异")
    if is_multiallelic.any():
        logger.warning("VCF 包含多等位基因位点")

    complex_mask = has_indel | is_multiallelic
    if complex_mask.any() and skip_complex_variants:
        logger.info("已启用 --skip-complex-variants: 仅保留二等位 SNP 进行计算")
        ref = vcf["variants/REF"].astype(str)
        alt = vcf["variants/ALT"]
        numalt = vcf["variants/numalt"]
        alt_len = np.vectorize(lambda x: len(x) if x is not None else 0)(alt)
        first_alt_len = alt_len[:, 0] if alt_len.ndim > 1 else alt_len
        is_snp = (np.char.str_len(ref) == 1) & (first_alt_len == 1)
        is_biallelic = numalt == 1
        chrom_mask = chrom_mask & is_snp & is_biallelic

    filtered: Dict[str, np.ndarray] = {}
    for key, value in vcf.items():
        if key == "samples":
            filtered[key] = value
        elif key.startswith("variants/"):
            filtered[key] = value[chrom_mask]
        elif key.startswith("calldata/"):
            filtered[key] = value[chrom_mask]
        else:
            if value.shape[0] == chrom_mask.shape[0]:
                filtered[key] = value[chrom_mask]
            else:
                filtered[key] = value

    if filtered["variants/POS"].size == 0:
        raise ValueError("过滤后无可用位点，请检查染色体或过滤条件")

    logger.info("过滤完成: 保留位点数 %s", filtered["variants/POS"].size)

    return filtered


def validate_sample_match(samples: np.ndarray, sample_df: pd.DataFrame, id_col: str) -> None:
    vcf_samples = set(samples.tolist())
    meta_samples = set(sample_df[id_col].astype(str).tolist())

    missing_in_vcf = sorted(meta_samples - vcf_samples)
    missing_in_meta = sorted(vcf_samples - meta_samples)

    logger.info("样本匹配校验结果")
    logger.info("VCF 样本数: %s", len(vcf_samples))
    logger.info("样本信息表样本数: %s", len(meta_samples))
    logger.info("样本信息表中未在 VCF 出现的样本: %s", len(missing_in_vcf))
    logger.info("VCF 中未在样本表出现的样本: %s", len(missing_in_meta))

    if missing_in_vcf or missing_in_meta:
        logger.warning("样本ID存在不匹配，将继续计算")


def build_populations(
    samples: np.ndarray,
    sample_df: pd.DataFrame,
    id_col: str,
    group_col: str,
) -> Dict[str, List[int]]:
    logger.info("构建群体索引: %s", group_col)
    sample_to_index = {sample: i for i, sample in enumerate(samples)}
    populations: Dict[str, List[int]] = {}

    for group_name, group_df in sample_df.groupby(group_col):
        if pd.isna(group_name):
            continue
        indices = []
        for sample_id in group_df[id_col].astype(str):
            if sample_id in sample_to_index:
                indices.append(sample_to_index[sample_id])
            else:
                logger.warning("样本 %s 在 VCF 中未找到", sample_id)

        if indices:
            populations[str(group_name)] = indices

    return populations
