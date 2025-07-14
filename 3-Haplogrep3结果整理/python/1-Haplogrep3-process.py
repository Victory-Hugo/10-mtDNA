import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from aquarel import load_theme
from itertools import chain
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LinearSegmentedColormap

def run_qc(
    input_tsv : str  = "../data/质量控制.tsv",
    output_dir: str  = "../output"
) -> None:
    """
    对 Haplogrep3 分型结果进行 mtDNA 质量控制、统计与可视化，所有输出写入 `output_dir`.

    参数
    ----
    input_tsv : str
        Haplogrep3 导出的 .tsv 文件路径。
    output_dir : str
        结果文件输出目录。

    主要输出
    --------
    ├─ 质量描述.csv
    ├─ quality.pdf
    ├─ 总突变描述.csv
    ├─ 变异计数综合图.pdf
    ├─ 非同义蛋白质分布.csv
    ├─ 非同义密码子分布.csv
    ├─ 非同义编码区突变分布.csv
    ├─ 非同义突变分布综合图.pdf
    ├─ 新发变异类型统计.csv
    ├─ 新发变异区域分布.csv
    ├─ 新发变异.pdf
    └─ 新发变异综合图.pdf
    """
    #--------------------------------------------------#
    # 0. 依赖检查 & 环境准备
    #--------------------------------------------------#
    os.makedirs(output_dir, exist_ok=True)
    theme = load_theme("boxy_light")

    #--------------------------------------------------#
    # 1. 预处理：替换 Remaining_Polys 中的 " (" → "("
    #--------------------------------------------------#
    df = pd.read_csv(input_tsv, sep="\t", encoding="utf-8", quotechar='"')
    df["Remaining_Polys"] = df["Remaining_Polys"].str.replace(" (", "(", regex=False)
    df.to_csv(input_tsv, sep="\t", index=False, encoding="utf-8", quotechar='"')

    #--------------------------------------------------#
    # 2. Quality 分布
    #--------------------------------------------------#
    df["Quality"].describe().to_csv(
        os.path.join(output_dir, "质量描述.csv"), encoding="utf-8"
    )
    theme.apply()
    plt.rcParams.update({"font.family": "Arial", "pdf.fonttype": 42, "ps.fonttype": 42})

    plt.figure(figsize=(5, 5))
    sns.histplot(df["Quality"], kde=True, bins=30, color="#0FB797",
                 edgecolor="black", alpha=1)
    plt.title("Distribution of quality scores of new data", fontsize=16)
    plt.xlabel("Quality score", fontsize=14)
    plt.ylabel("Frequency",     fontsize=14)
    plt.grid(False)
    plt.tight_layout()
    theme.apply_transforms()
    plt.savefig(os.path.join(output_dir, "quality.pdf"))
    plt.close()

    #--------------------------------------------------#
    # 3. 总突变计数（插入 / 缺失 / SNP）
    #--------------------------------------------------#
    all_mut = [m for subs in df["Input_Sample"].dropna().str.split()
                 for m in subs]

    insertions = [m for m in all_mut if "." in m]
    deletions  = [m for m in all_mut if "d" in m]
    snps       = [m for m in all_mut if "." not in m and "d" not in m]

    mut_df = pd.DataFrame({
        "Mutation type": ["Insertions", "Deletions", "SNVs"],
        "Total count"  : [len(insertions), len(deletions), len(snps)],
        "Unique count" : [len(set(insertions)), len(set(deletions)), len(set(snps))]
    })
    mut_df.to_csv(os.path.join(output_dir, "总突变描述.csv"), index=False)

    # 3-1. 总突变条形图
    melted = mut_df.melt(id_vars="Mutation type",
                         value_vars=["Total count", "Unique count"],
                         var_name="Count Type", value_name="Count")
    palette = ["#1A8986", "#0FB797", "#3DE0C0"]
    def _plot_tot(ax, data, title):
        sns.barplot(data=data, x="Mutation type", y="Count",
                    hue="Mutation type", palette=palette,
                    legend=False, ax=ax)
        ax.set_title(title, fontsize=14)
        ax.set_xlabel("Mutation type"); ax.set_ylabel("Count")
        ax.grid(False)
        for p in ax.patches:
            if p.get_height()>0:
                ax.annotate(f"{int(p.get_height())}",
                            (p.get_x()+p.get_width()/2, p.get_height()),
                            ha="center", va="bottom", xytext=(0,6),
                            textcoords="offset points", fontsize=9)

    fig, axes = plt.subplots(1, 2, figsize=(11.69, 11.69/2), constrained_layout=True)
    _plot_tot(axes[0], melted[melted["Count Type"]=="Total count"],  "Total mutation counts")
    _plot_tot(axes[1], melted[melted["Count Type"]=="Unique count"], "Unique mutation counts")
    theme.apply_transforms()
    with PdfPages(os.path.join(output_dir, "变异计数综合图.pdf")) as pdf:
        pdf.savefig(fig)
    plt.close(fig)

    #--------------------------------------------------#
    # 4. 非同义突变：蛋白质 / 密码子 / 编码区位置
    #--------------------------------------------------#
    df["AAC_In_Remainings"] = df["AAC_In_Remainings"].str.replace(" [", "[", regex=False)
    df["AAC_Split"] = df["AAC_In_Remainings"].dropna().apply(
        lambda x: [s.strip()+"]" for s in x.split("] ") if s]
    )
    aac_entries = list(chain.from_iterable(df["AAC_Split"].dropna()))
    aac_series  = pd.Series(aac_entries).value_counts()
    aac_df      = aac_series.reset_index()
    aac_df.columns = ["AAC_Entry", "Count"]

    # 辅助提取函数
    extract = lambda pat, x: (m.group(1) if (m:=re.search(pat,x)) else None)
    aac_df["Protein"]         = aac_df["AAC_Entry"].apply(lambda x: extract(r"\| ([^\|]+?) \]", x))
    aac_df["Codon_Position"]  = aac_df["AAC_Entry"].apply(lambda x: extract(r"\| ([^\|]+?)\|", x))
    aac_df["Numeric_Position"]= aac_df["AAC_Entry"].apply(lambda x: int(extract(r"^(\d+)", x) or 0))

    protein_df = aac_df.groupby("Protein")["Count"].sum().sort_values(ascending=False).reset_index()
    codon_df   = aac_df.groupby("Codon_Position")["Count"].sum().sort_values(ascending=False).reset_index()
    numpos_df  = aac_df.groupby("Numeric_Position")["Count"].sum().sort_values(ascending=False).reset_index()

    protein_df.to_csv(os.path.join(output_dir, "非同义蛋白质分布.csv"), index=False)
    codon_df  .to_csv(os.path.join(output_dir, "非同义密码子分布.csv"), index=False)
    numpos_df .to_csv(os.path.join(output_dir, "非同义编码区突变分布.csv"), index=False)

    # 4-1. 可视化
    codon_palette   = ["#DD6617", "#FF9300", "#FFB45A"]
    custom_colors   = ["#50184E","#B86265","#D49C87","#DCDFD2",
                       "#96B89B","#4E9280","#31646C","#003936"]
    cmap            = LinearSegmentedColormap.from_list("custom_gradient", custom_colors,
                                                        N=len(protein_df))
    protein_palette = [cmap(i/len(protein_df)) for i in range(len(protein_df))]

    def _plot_bar(ax, data, x, y, title, palette, rot=0):
        sns.barplot(data=data, x=x, y=y, hue=x, palette=palette,
                    legend=False, ax=ax)
        ax.set_title(title); ax.set_xlabel(x); ax.set_ylabel(y)
        ax.grid(False)
        ax.tick_params(axis="x", rotation=rot, labelsize=7)
        for p in ax.patches:
            if p.get_height()>0:
                ax.annotate(f"{int(p.get_height())}",
                            (p.get_x()+p.get_width()/2, p.get_height()),
                            ha="center", va="bottom", xytext=(0,3),
                            textcoords="offset points", fontsize=7)

    fig, axes = plt.subplots(1, 2, figsize=(11.69, 11.69/2), constrained_layout=True)
    _plot_bar(axes[0], codon_df,   "Codon_Position",  "Count",
              "Codon Position Distribution", codon_palette)
    _plot_bar(axes[1], protein_df, "Protein", "Count",
              "Proteins affected by non-synonymous mutations",
              protein_palette, rot=90)
    theme.apply_transforms()
    with PdfPages(os.path.join(output_dir, "非同义突变分布综合图.pdf")) as pdf:
        pdf.savefig(fig)
    plt.close(fig)

    #--------------------------------------------------#
    # 5. 新发变异 Remaining_Polys 统计
    #--------------------------------------------------#
    remain = df["Remaining_Polys"].dropna().str.split(expand=True).stack()
    total_vals, unique_vals = len(remain), remain.nunique()
    brackets   = remain.str.extract(r"\((.*?)\)")[0]
    value_pre  = remain.str.extract(r"(^[^\(]+)")[0]

    counts = {
        "Hotspot"              : brackets.eq("hotspot").sum(),
        "Unique hotspot"       : remain[brackets=="hotspot"].nunique(),
        "Global private"       : brackets.eq("globalPrivateMutation").sum(),
        "Unique global private": remain[brackets=="globalPrivateMutation"].nunique(),
        "Local private"        : brackets.eq("localPrivateMutation").sum(),
        "Unique local private" : remain[brackets=="localPrivateMutation"].nunique(),
        "Insertions"           : value_pre.str.contains(r"\.").sum(),
        "Deletions"            : value_pre.str.contains(r"d$").sum(),
    }
    counts["SNP"] = total_vals - counts["Insertions"] - counts["Deletions"]

    summary = pd.DataFrame({
        "Metric": ["Total values","Unique values"] + list(counts.keys()),
        "Count" : [total_vals, unique_vals] + list(counts.values())
    })
    summary.to_csv(os.path.join(output_dir, "新发变异类型统计.csv"), index=False)

    # 5-1. 区域分布
    def _bytype(mask):  # 提取位点并计数
        return value_pre[mask].str.extract(r"^(\d+)")[0].value_counts()
    region_df = pd.DataFrame({
        "Insertions": _bytype(value_pre.str.contains(r"\.")),
        "Deletions" : _bytype(value_pre.str.contains(r"d$")),
        "SNPs"      : _bytype(~value_pre.str.contains(r"\.|d$"))
    }).fillna(0).astype(int)
    region_df.to_csv(os.path.join(output_dir, "新发变异区域分布.csv"), index_label="Site")

    # 5-2. 新发变异条形图 (Unique / Total & Top 30 位置)
    unique_stats = summary[summary["Metric"].str.contains("Unique")]
    total_stats  = summary[~summary["Metric"].str.contains("Unique")]

    def _plot_stats(ax, data, title):
        sns.barplot(data=data, x="Metric", y="Count",
                    hue="Metric", palette=custom_colors[:len(data)],
                    legend=False, ax=ax)
        ax.set_title(title); ax.set_xlabel("Metric"); ax.set_ylabel("Count")
        ax.grid(False)
        for p in ax.patches:
            if p.get_height()>0:
                ax.annotate(f"{int(p.get_height())}",
                            (p.get_x()+p.get_width()/2, p.get_height()),
                            ha="center", va="bottom",
                            xytext=(0,3), textcoords="offset points", fontsize=8)

    fig, axes = plt.subplots(1, 2, figsize=(11.69, 11.69/2), constrained_layout=True)
    _plot_stats(axes[0], unique_stats, "Unique variant counts")
    _plot_stats(axes[1], total_stats,  "Total variant counts")
    theme.apply_transforms()
    with PdfPages(os.path.join(output_dir, "新发变异.pdf")) as pdf:
        pdf.savefig(fig)
    plt.close(fig)

    # Top30 位点
    reg = pd.read_csv(os.path.join(output_dir, "新发变异区域分布.csv"))
    top_ins = reg.sort_values("Insertions", ascending=False).head(30)
    top_del = reg.sort_values("Deletions",  ascending=False).head(30)
    top_snp = reg.sort_values("SNPs",       ascending=False).head(30)

    def _simple_bar(ax, df, col, title, color):
        bars = ax.bar(df["Site"].astype(str), df[col], color=color)
        ax.set_title(title); ax.set_xlabel("Site"); ax.set_ylabel(col)
        ax.tick_params(axis="x", rotation=90, labelsize=6)
        for b in bars:
            if b.get_height()>0:
                ax.text(b.get_x()+b.get_width()/2, b.get_height(),
                        f"{int(b.get_height())}", ha="center", va="bottom",
                        fontsize=5)

    fig, axes = plt.subplots(1, 3, figsize=(11.69, 11.69/3), constrained_layout=True)
    _simple_bar(axes[0], top_ins, "Insertions", "Top 30 sites by insertions", "#50184E")
    _simple_bar(axes[1], top_del, "Deletions",  "Top 30 sites by deletions",  "#003936")
    _simple_bar(axes[2], top_snp, "SNPs",       "Top 30 sites by SNPs",       "#00AFFF")
    theme.apply_transforms()
    with PdfPages(os.path.join(output_dir, "新发变异综合图.pdf")) as pdf:
        pdf.savefig(fig)
    plt.close(fig)

    print(f"全部分析完成，结果已保存至：{os.path.abspath(output_dir)}")

# ===================== 调用示例 =====================
run_qc("/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/3-Haplogrep3结果整理/data/20K.txt", "/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/3-Haplogrep3结果整理/output")