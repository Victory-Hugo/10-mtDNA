#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
import plotly.express as px

def main():
    parser = argparse.ArgumentParser(description="对 CSV 做 PCA 并输出结果（可指定降维维度）")
    parser.add_argument("--input", "-i", required=True,
                        help="输入 CSV 文件路径（第一列群体或样本标签，后续列为数值特征）")
    parser.add_argument("--output-dir", "-o", required=True,
                        help="输出目录，若不存在会自动创建")
    parser.add_argument("--n-components", "-n", type=int, default=20,
                        help="PCA 要保留的主成分个数（默认 20）")
    args = parser.parse_args()

    input_csv    = args.input
    out_dir      = args.output_dir
    n_components = args.n_components

    # 自动创建输出目录
    os.makedirs(out_dir, exist_ok=True)

    # 加载数据
    data = pd.read_csv(input_csv)
    # 第一列是群体或样本标签
    labels       = data.iloc[:, 0]
    numeric_data = data.iloc[:, 1:]

    # 执行 PCA
    pca = PCA(n_components=n_components)
    pcs = pca.fit_transform(numeric_data)
    evr = pca.explained_variance_ratio_

    # 1. 保存解释度
    evr_file = os.path.join(out_dir, f"explained_variance_{n_components}.csv")
    with open(evr_file, "w", encoding="utf-8") as f:
        f.write(f"Explained Variance Ratio ({n_components} components):\n")
        for i, r in enumerate(evr, start=1):
            f.write(f"PC{i},{r:.6f}\n")

    # 2. 保存 PCA 结果
    cols = [f"PC{i}" for i in range(1, n_components+1)]
    df_pca = pd.DataFrame(pcs, columns=cols)
    df_pca.insert(0, data.columns[0], labels)
    results_file = os.path.join(out_dir, f"pca_results_{n_components}.csv")
    df_pca.to_csv(results_file, sep=",", index=False, encoding="utf-8")

    # 3. 用 Seaborn 画前两主成分散点图
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x="PC1", y="PC2", hue=data.columns[0], data=df_pca)
    plt.title(f"PCA (PC1 vs PC2, {n_components} components)")
    plt.xlabel("PC1")
    plt.ylabel("PC2")
    plt.legend(title=data.columns[0], bbox_to_anchor=(1.05, 1), loc="upper left")
    plot_file = os.path.join(out_dir, f"pca_plot_{n_components}.png")
    plt.tight_layout()
    plt.savefig(plot_file, dpi=300)
    plt.close()

    print("PCA 完成，结果保存在：")
    print(f" - 解释度：{evr_file}")
    print(f" - PCA 坐标：{results_file}")
    print(f" - 散点图：{plot_file}")

if __name__ == "__main__":
    main()
