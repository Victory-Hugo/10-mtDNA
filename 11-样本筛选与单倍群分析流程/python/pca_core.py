#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""对频率矩阵执行 PCA 并输出解释度、坐标与散点图。"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Iterable

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from sklearn.decomposition import PCA

log = logging.getLogger(__name__)


def configure_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level.upper(), logging.INFO),
        format="[%(levelname)s] %(message)s",
    )


def run(input: str, output_dir: str, n_components: int = 20, log_level: str = "INFO") -> int:
    configure_logging(log_level)

    input_path = Path(input)
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    data = pd.read_csv(input_path)
    labels = data.iloc[:, 0]
    numeric_data = data.iloc[:, 1:]

    pca = PCA(n_components=n_components)
    pcs = pca.fit_transform(numeric_data)
    evr = pca.explained_variance_ratio_

    evr_file = out_dir / f"explained_variance_{n_components}.csv"
    with evr_file.open("w", encoding="utf-8") as handle:
        handle.write(f"Explained Variance Ratio ({n_components} components):\n")
        for idx, ratio in enumerate(evr, start=1):
            handle.write(f"PC{idx},{ratio:.6f}\n")

    cols = [f"PC{i}" for i in range(1, n_components + 1)]
    df_pca = pd.DataFrame(pcs, columns=cols)
    df_pca.insert(0, data.columns[0], labels)
    results_file = out_dir / f"pca_results_{n_components}.csv"
    df_pca.to_csv(results_file, sep=",", index=False, encoding="utf-8")

    plt.figure(figsize=(8, 6))
    sns.scatterplot(x="PC1", y="PC2", hue=data.columns[0], data=df_pca)
    plt.title(f"PCA (PC1 vs PC2, {n_components} components)")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend(title=data.columns[0], bbox_to_anchor=(1.05, 1), loc="upper left")
    plot_file = out_dir / f"pca_plot_{n_components}.png"
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300)
    plt.close()

    log.info("PCA 完成，解释度：%s", evr_file)
    log.info("PCA 坐标：%s", results_file)
    log.info("散点图：%s", plot_file)
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", "-i", required=True, help="输入 CSV 文件路径")
    parser.add_argument("--output-dir", "-o", required=True, help="输出目录")
    parser.add_argument("--n-components", "-n", type=int, default=20, help="PCA 主成分数量")
    parser.add_argument("--log-level", default="INFO", help="日志级别")
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    return run(**vars(args))


if __name__ == "__main__":
    raise SystemExit(main())
