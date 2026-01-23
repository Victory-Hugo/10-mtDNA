from __future__ import annotations

import argparse
import subprocess
from pathlib import Path


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Tajima's D 显著性工具封装")
    parser.add_argument("--tool", required=True, help="Tajima 显著性工具脚本路径")
    parser.add_argument("--vcf", required=True, help="输入 VCF")
    parser.add_argument("--sample-table", required=True, help="样本信息表")
    parser.add_argument("--group-col", required=True, help="分群列")
    parser.add_argument("--id-col", required=True, help="样本ID列")
    parser.add_argument("--chrom", required=True, help="染色体名")
    parser.add_argument("--out", required=True, help="输出路径")
    parser.add_argument("--n-replicates", type=int, required=True, help="模拟次数")
    parser.add_argument("--length", type=int, required=True, help="序列长度")
    parser.add_argument("--ne-min", type=float, required=True, help="Ne 下限")
    parser.add_argument("--ne-max", type=float, required=True, help="Ne 上限")
    parser.add_argument("--mu-min", type=float, required=True, help="突变率下限")
    parser.add_argument("--mu-max", type=float, required=True, help="突变率上限")
    parser.add_argument("--seed", type=int, default=None, help="随机种子")
    return parser.parse_args()


def main() -> None:
    args = _parse_args()
    tool = Path(args.tool)

    cmd = [
        "python",
        str(tool),
        "--vcf",
        args.vcf,
        "--sample-table",
        args.sample_table,
        "--group-col",
        args.group_col,
        "--id-col",
        args.id_col,
        "--chrom",
        args.chrom,
        "--out",
        args.out,
        "--n-replicates",
        str(args.n_replicates),
        "--length",
        str(args.length),
        "--ne-min",
        str(args.ne_min),
        "--ne-max",
        str(args.ne_max),
        "--mu-min",
        str(args.mu_min),
        "--mu-max",
        str(args.mu_max),
    ]

    if args.seed is not None:
        cmd.extend(["--seed", str(args.seed)])

    subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
