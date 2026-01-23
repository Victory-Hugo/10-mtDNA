from __future__ import annotations

import argparse
import logging
import subprocess
from pathlib import Path
from typing import List

import allel
import numpy as np
import pandas as pd

from preprocess import (
    build_populations,
    check_group_counts,
    filter_variants,
    read_sample_table,
    read_vcf,
    summarize_groups,
    validate_sample_match,
    validate_sample_table,
)
from metrics import (
    bootstrap_diversity,
    calculate_dxy,
    calculate_fst,
    calculate_pi,
    calculate_sfs,
    calculate_tajima_d,
    calculate_theta_w,
    permutation_test_fst_dxy,
    prepare_population_allele_counts,
)
from output_writer import save_pairwise_metrics, save_population_metrics

logger = logging.getLogger(__name__)


class _ColorFormatter(logging.Formatter):
    _colors = {
        logging.DEBUG: "\033[36m",
        logging.INFO: "\033[32m",
        logging.WARNING: "\033[33m",
        logging.ERROR: "\033[31m",
        logging.CRITICAL: "\033[41m\033[97m",
    }
    _reset = "\033[0m"

    def format(self, record: logging.LogRecord) -> str:
        message = super().format(record)
        color = self._colors.get(record.levelno, "")
        if not color:
            return message
        return f"{color}{message}{self._reset}"


def configure_logging(log_path: Path, level: str = "INFO") -> None:
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_level = logging.getLevelName(level.upper())

    file_handler = logging.FileHandler(log_path, encoding="utf-8")
    file_handler.setFormatter(
        logging.Formatter("%(asctime)s | %(levelname)s | %(name)s | %(message)s")
    )

    console_handler = logging.StreamHandler()
    console_handler.setFormatter(
        _ColorFormatter("%(asctime)s | %(levelname)s | %(message)s", datefmt="%H:%M:%S")
    )

    logging.basicConfig(level=log_level, handlers=[file_handler, console_handler], force=True)


def run(
    vcf_path: str,
    sample_table_path: str,
    group_cols: List[str],
    id_col: str,
    output_dir: str,
    chrom_name: str,
    skip_complex_variants: bool = False,
    enable_significance: bool = False,
    permutation_n: int = 1000,
    bootstrap_n: int = 1000,
    random_seed: int | None = None,
    tajima_tool: str | None = None,
    tajima_n_replicates: int = 2000,
    tajima_length: int = 16569,
    tajima_ne_min: float = 2000,
    tajima_ne_max: float = 20000,
    tajima_mu_min: float = 1e-8,
    tajima_mu_max: float = 3e-8,
) -> None:
    vcf_path = Path(vcf_path)
    sample_table_path = Path(sample_table_path)
    output_dir = Path(output_dir)

    logger.info("开始群体遗传学计算")
    logger.info("VCF: %s", vcf_path)
    logger.info("样本表: %s", sample_table_path)
    logger.info("分群列: %s", group_cols)
    logger.info("输出目录: %s", output_dir)
    logger.info("显著性检验: %s", "开启" if enable_significance else "关闭")

    logger.info("开始读取样本信息表")
    sample_df = read_sample_table(sample_table_path)
    validate_sample_table(sample_df, id_col, group_cols)
    logger.info("样本信息表读取完成，样本数: %s", len(sample_df))

    logger.info("开始读取 VCF")
    vcf = read_vcf(vcf_path)
    logger.info("VCF 格式校验通过")

    validate_sample_match(vcf["samples"], sample_df, id_col)

    vcf = filter_variants(vcf, chrom_name, skip_complex_variants)
    logger.info("染色体 %s 校验通过，位点数: %s", chrom_name, len(vcf["variants/POS"]))

    gt = allel.GenotypeArray(vcf["calldata/GT"])
    pos = vcf["variants/POS"]

    vcf_samples = set(vcf["samples"].astype(str).tolist())
    sample_df = sample_df[sample_df[id_col].astype(str).isin(vcf_samples)].copy()
    logger.info("按 VCF 过滤样本表后样本数: %s", len(sample_df))

    for group_col in group_cols:
        logger.info("开始处理分群列: %s", group_col)
        group_counts = summarize_groups(sample_df, group_col)
        check_group_counts(group_counts, group_col)

        populations = build_populations(vcf["samples"], sample_df, id_col, group_col)
        if len(populations) < 2:
            raise ValueError(f"分群列 {group_col} 的有效群体数量不足 2，无法计算 FST/Dxy")

        logger.info("分群列 %s: 群体数 %s，开始计算指标", group_col, len(populations))

        logger.info("计算 π (nucleotide diversity)")
        pi_results = calculate_pi(gt, populations, pos)
        logger.info("计算 θw (Watterson's theta)")
        theta_w_results = calculate_theta_w(gt, populations, pos)
        logger.info("计算 Tajima's D")
        tajima_d_results = calculate_tajima_d(gt, populations, pos)
        logger.info("计算 folded SFS")
        sfs_folded_results = calculate_sfs(gt, populations, folded=True)
        logger.info("计算 unfolded SFS")
        sfs_unfolded_results = calculate_sfs(gt, populations, folded=False)
        logger.info("计算 FST")
        fst_results = calculate_fst(gt, populations)
        logger.info("计算 Dxy")
        dxy_results = calculate_dxy(gt, populations, pos)

        fst_pvalues = None
        dxy_pvalues = None
        pi_bootstrap = None
        theta_w_bootstrap = None

        group_out = output_dir / f"group_{group_col}"

        if enable_significance:
            logger.info("显著性检验: 置换 FST/Dxy (n=%s)", permutation_n)
            fst_pvalues, dxy_pvalues = permutation_test_fst_dxy(
                gt,
                populations,
                pos,
                n_perm=permutation_n,
                seed=random_seed,
            )

            logger.info("显著性检验: bootstrap π/θw (n=%s)", bootstrap_n)
            ac_by_pop = prepare_population_allele_counts(gt, populations)
            boot = bootstrap_diversity(ac_by_pop, pos, n_boot=bootstrap_n, seed=random_seed)

            pi_bootstrap = _pairwise_bootstrap_table(
                boot,
                pi_results,
                metric_name="pi",
            )
            theta_w_bootstrap = _pairwise_bootstrap_table(
                boot,
                theta_w_results,
                metric_name="theta_w",
            )

            if tajima_tool:
                logger.info("显著性检验: Tajima's D 外部工具")
                _run_tajima_tool(
                    tajima_tool,
                    group_out / "population",
                    vcf_path,
                    sample_table_path,
                    group_col,
                    id_col,
                    chrom_name,
                    tajima_n_replicates,
                    tajima_length,
                    tajima_ne_min,
                    tajima_ne_max,
                    tajima_mu_min,
                    tajima_mu_max,
                    random_seed,
                )
            else:
                logger.warning("未设置 Tajima's D 外部工具，跳过显著性检验")

        save_population_metrics(
            group_out / "population",
            pi_results,
            theta_w_results,
            tajima_d_results,
            sfs_folded_results,
            sfs_unfolded_results,
            group_counts,
            pi_bootstrap=pi_bootstrap,
            theta_w_bootstrap=theta_w_bootstrap,
        )
        save_pairwise_metrics(
            group_out / "pairwise",
            fst_results,
            dxy_results,
            fst_p=fst_pvalues,
            dxy_p=dxy_pvalues,
        )

        logger.info("分群列 %s 计算完成，结果已写入: %s", group_col, group_out)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="群体遗传学计算模块 (scikit-allel)")
    parser.add_argument("--vcf", required=True, help="输入的 VCF(.vcf/.vcf.gz) 文件路径")
    parser.add_argument("--sample-table", required=True, help="样本信息表格路径 (.csv/.tsv)")
    parser.add_argument("--group-cols", required=True, help="分群列名，多个列用英文逗号分隔")
    parser.add_argument("--id-col", default="ID", help="样本ID列名，默认 ID")
    parser.add_argument("--output-dir", required=True, help="结果输出目录")
    parser.add_argument("--chrom", required=True, help="计算的染色体名称，如 chrM")
    parser.add_argument("--skip-complex-variants", action="store_true", help="跳过 InDel/多等位基因位点")
    parser.add_argument("--enable-significance", action="store_true", help="开启统计显著性检验")
    parser.add_argument("--permutation-n", type=int, default=1000, help="FST/Dxy 置换次数")
    parser.add_argument("--bootstrap-n", type=int, default=1000, help="π/θw bootstrap 次数")
    parser.add_argument("--random-seed", type=int, default=None, help="随机种子")
    parser.add_argument(
        "--tajima-tool",
        default=None,
        help="Tajima's D 显著性外部工具脚本路径",
    )
    parser.add_argument("--tajima-n-replicates", type=int, default=2000, help="Tajima 模拟次数")
    parser.add_argument("--tajima-length", type=int, default=16569, help="Tajima 序列长度")
    parser.add_argument("--tajima-ne-min", type=float, default=2000, help="Tajima Ne 下限")
    parser.add_argument("--tajima-ne-max", type=float, default=20000, help="Tajima Ne 上限")
    parser.add_argument("--tajima-mu-min", type=float, default=1e-8, help="Tajima 突变率下限")
    parser.add_argument("--tajima-mu-max", type=float, default=3e-8, help="Tajima 突变率上限")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    group_cols = [col.strip() for col in args.group_cols.split(",") if col.strip()]

    configure_logging(Path(args.output_dir) / "population_genetics.log")

    run(
        vcf_path=args.vcf,
        sample_table_path=args.sample_table,
        group_cols=group_cols,
        id_col=args.id_col,
        output_dir=args.output_dir,
        chrom_name=args.chrom,
        skip_complex_variants=args.skip_complex_variants,
        enable_significance=args.enable_significance,
        permutation_n=args.permutation_n,
        bootstrap_n=args.bootstrap_n,
        random_seed=args.random_seed,
        tajima_tool=args.tajima_tool,
        tajima_n_replicates=args.tajima_n_replicates,
        tajima_length=args.tajima_length,
        tajima_ne_min=args.tajima_ne_min,
        tajima_ne_max=args.tajima_ne_max,
        tajima_mu_min=args.tajima_mu_min,
        tajima_mu_max=args.tajima_mu_max,
    )


def _pairwise_bootstrap_table(
    boot: dict,
    observed: dict,
    metric_name: str,
) -> pd.DataFrame:
    pop_names = list(observed.keys())
    rows = []
    for i, pop1 in enumerate(pop_names):
        for j, pop2 in enumerate(pop_names):
            if i < j:
                obs_diff = float(observed[pop1] - observed[pop2])
                dist = boot[pop1][metric_name] - boot[pop2][metric_name]
                p_val = (1.0 + np.sum(np.abs(dist) >= abs(obs_diff))) / (len(dist) + 1.0)
                ci_low, ci_high = np.quantile(dist, [0.025, 0.975])
                rows.append(
                    {
                        "pop1": pop1,
                        "pop2": pop2,
                        "metric": metric_name,
                        "diff": obs_diff,
                        "ci_low": float(ci_low),
                        "ci_high": float(ci_high),
                        "p_value": float(p_val),
                    }
                )
    return pd.DataFrame(rows)


def _run_tajima_tool(
    tool_path: str,
    output_dir: Path,
    vcf_path: Path,
    sample_table_path: Path,
    group_col: str,
    id_col: str,
    chrom_name: str,
    n_replicates: int,
    length: int,
    ne_min: float,
    ne_max: float,
    mu_min: float,
    mu_max: float,
    seed: int | None,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "python",
        tool_path,
        "--vcf",
        str(vcf_path),
        "--sample-table",
        str(sample_table_path),
        "--group-col",
        group_col,
        "--id-col",
        id_col,
        "--chrom",
        chrom_name,
        "--out",
        str(output_dir / "tajima_significance.csv"),
        "--n-replicates",
        str(n_replicates),
        "--length",
        str(length),
        "--ne-min",
        str(ne_min),
        "--ne-max",
        str(ne_max),
        "--mu-min",
        str(mu_min),
        "--mu-max",
        str(mu_max),
    ]
    if seed is not None:
        cmd.extend(["--seed", str(seed)])
    logger.info("运行 Tajima 工具命令: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
