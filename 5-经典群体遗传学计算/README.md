# 经典群体遗传学计算（scikit-allel）

该模块读取 VCF 与样本信息表，按分群列计算：π、θw、Tajima's D、SFS、FST、Dxy。

## 目录结构

```
conf/
input/
data/
output/
tmp/
log/
script/
python/
requirements.txt
```

## 使用方式

```bash
python python/population_genetics.py \
  --vcf /path/to/input.vcf.gz \
  --sample-table /path/to/meta.tsv \
  --group-cols Population,Region \
  --id-col ID \
  --output-dir /path/to/output \
  --chrom chrM \
  --skip-complex-variants \
  --enable-significance \
  --permutation-n 1000 \
  --bootstrap-n 1000 \
  --random-seed 2026
```

## 输出说明

每个分群列会生成一个子目录：

```
output/
  group_Population/
    population/
      population_counts.csv
      pi.csv
      theta_w.csv
      tajima_d.csv
      sfs_long.csv
    pairwise/
      fst.csv
      dxy.csv
      fst_pvalue.csv
      dxy_pvalue.csv
      pi_bootstrap.csv
      theta_w_bootstrap.csv
```

## 指标公式核对（metrics.py）

以下核对基于 `python/metrics.py` 中实现，逐条对照 scikit-allel 官方函数定义：

- **π (nucleotide diversity)**
  - 实现：`allel.sequence_diversity(pos, ac)`
  - 含义：基于等位基因计数的每位点平均差异数，按位点坐标 `pos` 计算。
  - 结论：与 scikit-allel 定义一致，无明显逻辑偏差。

- **θw (Watterson's theta)**
  - 实现：`allel.watterson_theta(pos, ac)`
  - 含义：基于分离位点数与样本量的 $theta_w$，按位点坐标 `pos` 计算。
  - 结论：与 scikit-allel 定义一致，无明显逻辑偏差。

- **Tajima's D**
  - 实现：`allel.tajima_d(ac, pos=pos)`
  - 含义：由 $pi$ 与 $theta_w$ 的差异构建的检验统计量。
  - 结论：与 scikit-allel 定义一致；结果依赖于位点过滤和样本量。

- **SFS (site frequency spectrum)**
  - Folded：`allel.sfs_folded(ac_seg)`
  - Unfolded：`allel.sfs(dac)`，其中 `dac = ac_seg[:, 1]`
  - 结论：folded SFS 与 scikit-allel 一致。unfolded SFS 假设等位基因 1 为“导出等位”，若未指定外群或已定向，则属于“伪 unfold”。

- **FST (Weir & Cockerham)**
  - 实现：`allel.weir_cockerham_fst(gt_array, subpops)`，最终值为 $\sum a / \sum(a+b+c)$ 并截断到 $[0,1]$。
  - 结论：与 scikit-allel 推荐汇总方式一致。

- **Dxy (sequence divergence)**
  - 实现：`allel.sequence_divergence(pos, ac1, ac2)`
  - 含义：两群体间每位点的平均序列差异。
  - 结论：与 scikit-allel 定义一致。

- **Bootstrap (π/θw)**
  - 实现：对位点进行有放回重采样，分别调用 `sequence_diversity` 与 `watterson_theta`。
  - 结论：方法符合常规位点 bootstrap。

- **Permutation (FST/Dxy)**
  - 实现：对两群体样本标签置换，重复计算 FST 与 Dxy。
  - 结论：逻辑正确；适用于检验组间差异显著性。

**已知假设与注意**
- 结果依赖于 `filter_variants` 的位点过滤（如仅二等位 SNP）。
- 缺失基因型通过 `count_alleles()` 自动忽略；若缺失较多，建议报告有效样本量。
- Unfolded SFS 的“导出等位”需要外群定向，否则只能解读为相对频谱。

## Tajima's D 显著性（msprime）

推荐使用 `msprime` 模拟中性模型生成 Tajima's D 的零分布，输出 p 值：

推荐在 `conf/1-run.conf` 中用参数化方式配置，避免长命令：

```bash
TAJIMA_TOOL="${PROJECT_DIR}/python/tajima_significance.py"
TAJIMA_N_REPLICATES="2000"
TAJIMA_LENGTH="16569"
TAJIMA_NE_MIN="2000"
TAJIMA_NE_MAX="20000"
TAJIMA_MU_MIN="1e-8"
TAJIMA_MU_MAX="3e-8"
```

输出字段包含 `population`, `tajima_d`, `p_value`, `n_samples`, `n_segregating` 以及模拟参数范围与 null CI。
```
