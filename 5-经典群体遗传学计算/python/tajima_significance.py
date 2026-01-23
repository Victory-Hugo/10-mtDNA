from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Dict, List, Tuple

import allel
import msprime
import numpy as np
import pandas as pd

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


def configure_logging(level: str = "INFO") -> None:
    log_level = logging.getLevelName(level.upper())
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(
        _ColorFormatter("%(asctime)s | %(levelname)s | %(message)s", datefmt="%H:%M:%S")
    )
    logging.basicConfig(level=log_level, handlers=[console_handler], force=True)


def read_sample_table(table_path: Path) -> pd.DataFrame:
    if table_path.suffix.lower() == ".tsv":
        return pd.read_csv(table_path, sep="\t")
    if table_path.suffix.lower() == ".csv":
        return pd.read_csv(table_path)
    raise ValueError("样本信息表格仅支持 .csv 或 .tsv 格式")


def validate_sample_table(sample_df: pd.DataFrame, id_col: str, group_col: str) -> None:
    missing_cols = [col for col in [id_col, group_col] if col not in sample_df.columns]
    if missing_cols:
        raise ValueError(f"样本信息表格缺少列: {missing_cols}")


def build_populations(
    samples: np.ndarray,
    sample_df: pd.DataFrame,
    id_col: str,
    group_col: str,
) -> Dict[str, List[int]]:
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


def read_vcf(vcf_path: Path, chrom: str | None) -> Tuple[np.ndarray, np.ndarray, allel.GenotypeArray]:
    vcf = allel.read_vcf(
        str(vcf_path),
        fields=["samples", "variants/CHROM", "variants/POS", "calldata/GT"],
    )
    samples = vcf["samples"]
    chroms = vcf["variants/CHROM"].astype(str)
    pos = vcf["variants/POS"]
    gt = allel.GenotypeArray(vcf["calldata/GT"])

    if chrom:
        mask = chroms == chrom
        pos = pos[mask]
        gt = gt.compress(mask, axis=0)

    if pos.size == 0:
        raise ValueError("过滤后无可用位点，请检查染色体")

    return samples, pos, gt


def observed_tajima_d(
    gt: allel.GenotypeArray,
    pos: np.ndarray,
    populations: Dict[str, List[int]],
) -> Dict[str, Dict[str, float]]:
    results = {}
    for pop, indices in populations.items():
        pop_gt = gt.take(indices, axis=1)
        ac = pop_gt.count_alleles()
        d_value = allel.tajima_d(ac, pos=pos)
        segregating = np.sum(np.sum(ac > 0, axis=1) > 1)
        results[pop] = {
            "tajima_d": float(d_value) if np.isfinite(d_value) else np.nan,
            "n_samples": int(pop_gt.n_samples),
            "n_segregating": int(segregating),
        }
    return results


def simulate_tajima_d(
    n_samples: int,
    n_reps: int,
    sequence_length: int,
    ne_min: float,
    ne_max: float,
    mu_min: float,
    mu_max: float,
    seed: int | None,
) -> np.ndarray:
    rng = np.random.default_rng(seed)
    results = []
    attempts = 0
    max_attempts = n_reps * 10

    while len(results) < n_reps and attempts < max_attempts:
        ne = float(rng.uniform(ne_min, ne_max))
        mu = float(rng.uniform(mu_min, mu_max))
        seed_anc = int(rng.integers(1, 2**31 - 1))
        seed_mut = int(rng.integers(1, 2**31 - 1))

        ts = msprime.sim_ancestry(
            samples=n_samples,
            ploidy=1,
            sequence_length=sequence_length,
            recombination_rate=0.0,
            population_size=ne,
            random_seed=seed_anc,
        )
        mts = msprime.sim_mutations(ts, rate=mu, random_seed=seed_mut)

        if mts.num_sites == 0:
            attempts += 1
            continue

        gm = mts.genotype_matrix()
        if gm.size == 0:
            attempts += 1
            continue

        gt = allel.GenotypeArray(gm[:, :, None])
        ac = gt.count_alleles()
        pos = mts.tables.sites.position
        d_value = allel.tajima_d(ac, pos=pos)
        if np.isfinite(d_value):
            results.append(float(d_value))
        attempts += 1

    return np.asarray(results, dtype=float)


def build_result_rows(
    obs: Dict[str, Dict[str, float]],
    sim: Dict[str, np.ndarray],
    sequence_length: int,
    ne_min: float,
    ne_max: float,
    mu_min: float,
    mu_max: float,
    n_reps: int,
) -> pd.DataFrame:
    rows = []
    for pop, info in obs.items():
        sim_values = sim.get(pop, np.array([]))
        if sim_values.size == 0 or not np.isfinite(info["tajima_d"]):
            rows.append(
                {
                    "population": pop,
                    "tajima_d": info["tajima_d"],
                    "p_value": np.nan,
                    "n_replicates": int(sim_values.size),
                    "n_samples": info["n_samples"],
                    "n_segregating": info["n_segregating"],
                    "null_ci_low": np.nan,
                    "null_ci_high": np.nan,
                    "length": sequence_length,
                    "ne_min": ne_min,
                    "ne_max": ne_max,
                    "mu_min": mu_min,
                    "mu_max": mu_max,
                }
            )
            continue

        p_value = (1.0 + np.sum(np.abs(sim_values) >= abs(info["tajima_d"]))) / (len(sim_values) + 1.0)
        ci_low, ci_high = np.quantile(sim_values, [0.025, 0.975])

        rows.append(
            {
                "population": pop,
                "tajima_d": info["tajima_d"],
                "p_value": float(p_value),
                "n_replicates": int(sim_values.size),
                "n_samples": info["n_samples"],
                "n_segregating": info["n_segregating"],
                "null_ci_low": float(ci_low),
                "null_ci_high": float(ci_high),
                "length": sequence_length,
                "ne_min": ne_min,
                "ne_max": ne_max,
                "mu_min": mu_min,
                "mu_max": mu_max,
            }
        )

    return pd.DataFrame(rows)


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Tajima's D 显著性检验 (msprime 模拟)")
    parser.add_argument("--vcf", required=True, help="输入 VCF 文件")
    parser.add_argument("--sample-table", required=True, help="样本信息表")
    parser.add_argument("--group-col", required=True, help="分群列名")
    parser.add_argument("--id-col", default="ID", help="样本ID列名")
    parser.add_argument("--chrom", default=None, help="染色体名，如 chrM")
    parser.add_argument("--out", required=True, help="输出 CSV 路径")
    parser.add_argument("--n-replicates", type=int, default=2000, help="模拟次数")
    parser.add_argument("--length", type=int, default=16569, help="序列长度 (bp)")
    parser.add_argument("--ne-min", type=float, default=2000, help="Ne 下限")
    parser.add_argument("--ne-max", type=float, default=20000, help="Ne 上限")
    parser.add_argument("--mu-min", type=float, default=1e-8, help="突变率下限 (per site per generation)")
    parser.add_argument("--mu-max", type=float, default=3e-8, help="突变率上限 (per site per generation)")
    parser.add_argument("--seed", type=int, default=2026, help="随机种子")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    configure_logging()

    vcf_path = Path(args.vcf)
    sample_table_path = Path(args.sample_table)
    output_path = Path(args.out)

    logger.info("读取 VCF 与样本表")
    sample_df = read_sample_table(sample_table_path)
    validate_sample_table(sample_df, args.id_col, args.group_col)
    samples, pos, gt = read_vcf(vcf_path, args.chrom)

    populations = build_populations(samples, sample_df, args.id_col, args.group_col)
    if not populations:
        raise ValueError("未找到可用群体")

    logger.info("计算观测 Tajima's D")
    obs = observed_tajima_d(gt, pos, populations)

    logger.info("模拟中性模型分布")
    sim_values: Dict[str, np.ndarray] = {}
    for pop, info in obs.items():
        n_samples = info["n_samples"]
        logger.info("群体 %s: n_samples=%s", pop, n_samples)
        sim_values[pop] = simulate_tajima_d(
            n_samples=n_samples,
            n_reps=args.n_replicates,
            sequence_length=args.length,
            ne_min=args.ne_min,
            ne_max=args.ne_max,
            mu_min=args.mu_min,
            mu_max=args.mu_max,
            seed=args.seed,
        )

    result = build_result_rows(
        obs=obs,
        sim=sim_values,
        sequence_length=args.length,
        ne_min=args.ne_min,
        ne_max=args.ne_max,
        mu_min=args.mu_min,
        mu_max=args.mu_max,
        n_reps=args.n_replicates,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(output_path, index=False)
    logger.info("结果已写入: %s", output_path)


if __name__ == "__main__":
    main()
