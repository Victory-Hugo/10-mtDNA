# scikit-allel Diversity Pipeline

该项目按 `pipeline-coding-standard` 组织，使用 `pipe/pipeline.sh` 调度多个双模式 Python 模块。

当前支持两类运行模式：
- 常规统计：默认启用
- 统一降采样 + 样本 bootstrap：由 `conf/Config.json` 中 `resampling.enable` 显式控制，默认关闭

## 运行

```bash
cd /mnt/c/Users/Administrator/Desktop/scikit-allel
conda create -y -n ScikitAllele -c conda-forge --override-channels python=3.10 numpy=1.26 pandas pyyaml zarr=2.18 numcodecs pip
conda run -n ScikitAllele pip install -r requirements.txt
bash pipe/pipeline.sh conf/Config.json
```

若要运行 demo bootstrap：

```bash
bash pipe/pipeline.sh conf/Config.bootstrap_demo.json
```

## 输入

- `data/merged_clean.vcf.gz`
- `meta/Population_Count20.list.tsv`

## ⚙️ 配置说明

主配置文件是 `conf/Config.json`。`conf/Config.bootstrap_demo.json` 是一个开启 bootstrap 的示例配置，具体重复次数以其中的 `bootstrap_replicates` 为准。

### 🧾 顶层字段

- `project_name`
  用于命名中间 Zarr 目录等产物。
- `input_vcf`
  输入 VCF 文件路径，建议使用 `bgzip` 压缩的 `.vcf.gz`。
- `sample_table`
  样本分组表路径，至少需要样本 ID 列和群体列。
- `sample_id_column`
  样本 ID 列名，默认示例为 `ID`。
- `group_column`
  群体列名，默认示例为 `Group`。
- `min_group_size`
  预处理时保留群体的最小样本数，小于该值的群体会被过滤。
- `contig_name`
  要分析的 contig 名称。当前示例是 `chrM`。
- `sequence_length`
  分析区域长度。对 `pi`、`theta_w`、`dxy` 这类按长度归一化的指标很重要。

### 📏 `sequence_length` 是否必须？

- 从“程序能否运行”的角度：不是绝对必须。
- 从“统计结果是否合理”的角度：强烈建议填写真实长度。
- 当前实现中，如果 `sequence_length <= 0`，程序会退回到“该 contig 的最大变异位点坐标”作为长度。
- 这会影响多样性指标的分母，可能造成系统偏差。
- 对 mtDNA 示例，建议继续使用 `16569`。

### 📂 路径相关字段

- `output_dir`
  输出目录根路径。
- `tmp_dir`
  中间文件目录，主要用于 Zarr。
- `meta_dir`
  元数据与预处理表格输出目录。
- `log_dir`
  日志目录。

### 🧬 `metrics`

- `within_groups`
  群体内指标列表。当前支持：
  - `pi`
  - `theta_w`
  - `tajima_d`
  - `haplotype_diversity`
- `between_groups`
  群体间指标列表。当前支持：
  - `dxy`
  - `hudson_fst`

### 🪟 `windows`

- `enable`
  是否开启滑窗统计，默认 `false`。
- `size`
  滑窗大小。
- `step`
  滑窗步长。

### 🔁 `resampling`

- `enable`
  是否开启统一降采样 + 样本 bootstrap。默认必须为 `false`。
- `sample_size_strategy`
  当前支持两种模式：
  - `min_group_size`：每个群体统一降采样到当前保留群体中的最小样本数。
  - `resampling_downsample`：每个群体固定按 `resampling_downsample` 指定的样本数进行有放回抽样。
- `resampling_downsample`
  当 `sample_size_strategy = resampling_downsample` 时生效。表示每次 bootstrap 时每个群体固定抽取多少个样本，默认可设为 `20`。
- `bootstrap_replicates`
  bootstrap 重复次数。常见正式分析可设为 `1000`，也可以按算力和输出体量调整。
- `random_seed`
  随机种子，用于保证 bootstrap 可复现。
- `write_replicate_tables`
  是否额外输出每一次 bootstrap 的明细结果。默认 `false`。开启后会生成 replicate-level TSV，便于查看“某个群体做了第几次抽样、这一次算出的数值是多少”。

### 🚀 `runtime`

- `python_bin`
  Python 解释器路径，当前指向 `ScikitAllele` conda 环境。
- `n_threads`
  并行线程或进程数。常规统计主要影响 bootstrap 并行。
- `overwrite_zarr`
  是否覆盖已有 Zarr 中间文件。
- `zarr_chunk_length`
  Zarr 按变异位点切块的长度。
- `zarr_chunk_width`
  Zarr 按样本方向切块的宽度。
- `alt_number`
  导入 VCF 时允许的 ALT 等位基因上限。

### ✅ 推荐用法

- 常规汇总分析：使用 `conf/Config.json`
- 示例 bootstrap 验证：使用 `conf/Config.bootstrap_demo.json`
- 正式 bootstrap 分析：把 `conf/Config.json` 中 `resampling.enable` 改为 `true`，再把 `bootstrap_replicates` 调到目标值

## 输出

- `output/0-qc/input_summary.tsv`
- `output/1-within/within_group_summary.tsv`
- `output/2-between/between_group_summary.tsv`
- `output/3-windows/windowed_metrics.tsv`
- `output/4-bootstrap/within_group_bootstrap_summary.tsv`
- `output/4-bootstrap/between_group_bootstrap_summary.tsv`
- `output/4-bootstrap/bootstrap_run_summary.tsv`
- `output/4-bootstrap/within_group_bootstrap_replicates.tsv`（仅当 `resampling.write_replicate_tables=true` 时生成）
- `output/4-bootstrap/between_group_bootstrap_replicates.tsv`（仅当 `resampling.write_replicate_tables=true` 时生成）

默认主配置关闭滑窗统计，也关闭 bootstrap；如需开启，可编辑 `conf/Config.json` 中的 `windows.enable` 和 `resampling.enable`。`conf/Config.bootstrap_demo.json` 提供了一个开启 bootstrap 的示例配置，是否额外输出每次 replicate 明细则由 `resampling.write_replicate_tables` 控制。
