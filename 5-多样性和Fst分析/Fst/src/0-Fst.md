# Fst计算脚本使用指南（增强版）

## 概述
这个增强版的Fst计算脚本基于PopGenome包，不仅计算配对Fst值，还提供**permutation检验**进行统计显著性检验，包括p值和置信区间。

## 功能特点
- 计算群体间配对Fst值
- **Permutation检验**进行统计显著性检验（更适合Fst统计量）
- 输出p值和基于permutation的置信区间
- 支持多种Fst计算方法
- 多核并行计算支持
- 详细的结果输出文件

## 统计方法说明

### 为什么使用Permutation而不是Bootstrap？
- **Permutation test**: 通过随机重排样本的群体标签来检验H0: Fst = 0，更适合检验群体分化的显著性
- **Bootstrap**: 重采样现有数据来估计参数的抽样分布，更适合估计置信区间
- 对于Fst统计量，permutation test是更合适的统计检验方法

### Permutation检验原理
1. 计算观察到的Fst值
2. 随机重排样本的群体标签N次
3. 每次重排后重新计算Fst值
4. p值 = (≥观察值的permutation次数 + 1) / (总permutation次数 + 1)

## 使用方法

### 基本语法
```bash
Rscript 0-Fst.R <fasta_file> <population_file> <output_file> [method] [cores] [permutations] [alpha]
```

### 参数说明
- `fasta_file`: 对齐的FASTA序列文件路径
- `population_file`: 群体信息文件路径（第一列：样本ID，第二列：群体名）
- `output_file`: 输出CSV文件路径
- `method`: 可选，Fst计算方法（默认：hudson）
- `cores`: 可选，使用的CPU核心数（默认：1）
- `permutations`: 可选，permutation检验次数（默认：1000，设为0跳过统计检验）
- `alpha`: 可选，显著性水平（默认：0.05）

### 支持的Fst计算方法
- `hudson`: Hudson et al. (1992) F_ST（默认）
- `nei`: Nei's G_ST
- `nucleotide`: 基于核苷酸的F_ST
- `haplotype`: 基于单倍型的F_ST

## 输出文件

### 基本输出
- `<output_file>`: 配对Fst矩阵（CSV格式）

### 统计检验输出（当permutations > 0时）
- `<output_file>_detailed_stats.csv`: 详细的配对比较表格，包含Fst值、p值、置信区间和显著性
- `<output_file>_pvalues.csv`: p值矩阵
- `<output_file>_confidence_intervals.csv`: 基于permutation的置信区间矩阵

## 使用示例

### 1. 基本Fst计算（无统计检验）
```bash
Rscript src/0-Fst.R example/Example.aln.fasta example/Example.txt output/fst_basic.csv hudson 1 0
```

### 2. 带统计检验的Fst计算（推荐）
```bash
Rscript src/0-Fst.R example/Example.aln.fasta example/Example.txt output/fst_with_stats.csv hudson 4 1000 0.05
```

### 3. 快速测试（少量permutation）
```bash
Rscript src/0-Fst.R example/Example.aln.fasta example/Example.txt output/fst_test.csv hudson 1 100 0.05
```

## 文件格式要求

### 群体信息文件格式
制表符分隔的文本文件，包含两列：
```
Sample_ID1    Population_A
Sample_ID2    Population_A
Sample_ID3    Population_B
Sample_ID4    Population_B
```

### FASTA文件格式
对齐的FASTA序列文件，样本ID需与群体文件中的ID匹配：
```
>Sample_ID1
ATCGATCGATCG...
>Sample_ID2
ATCGATCGATCG...
```

## 结果解释

### 详细统计文件字段说明
- `Population1`, `Population2`: 比较的两个群体
- `Fst`: 配对Fst值
- `P_value`: 统计显著性p值（H0: Fst = 0）
- `CI_lower`, `CI_upper`: 置信区间下限和上限
- `Significant`: 是否在给定显著性水平下显著

### 统计解释
- **Fst值**: 范围0-1，值越大表示群体分化越强
- **p值**: 小于alpha值（如0.05）表示统计显著
- **置信区间**: 包含真实Fst值的可能范围

## 注意事项

1. **计算时间**: Bootstrap检验会显著增加计算时间，特别是当样本数和bootstrap次数较大时
2. **内存使用**: 大型数据集可能需要较多内存
3. **数据质量**: 确保FASTA序列已正确对齐
4. **样本匹配**: 群体文件中的样本ID必须与FASTA文件中的序列ID完全匹配

## 推荐设置

### 快速测试
- Bootstrap: 100-500次
- 核心数: 1-2个

### 正式分析
- Bootstrap: 1000-10000次
- 核心数: 根据可用CPU调整
- 显著性水平: 0.05或0.01

## 故障排除

### 常见错误
1. **文件不存在**: 检查文件路径是否正确
2. **样本ID不匹配**: 确保群体文件和FASTA文件中的ID一致
3. **内存不足**: 减少bootstrap次数或使用更大内存的机器
4. **R包缺失**: 确保安装了PopGenome、parallel和tools包

### 安装依赖包
```r
install.packages(c("PopGenome", "parallel", "tools"))
```