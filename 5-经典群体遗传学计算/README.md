# 经典群体遗传学计算

该模块通过统一脚本 `pipe/1-run.sh` 读取配置文件，完成样本表过滤与主流程分析。脚本支持断点续跑、日志追踪、可选复杂变异跳过与显著性分析。

## 目录结构

```
conf/
log/
output/
python/
script/
```

## 快速开始

1. 编辑配置文件（示例：`conf/1-run.conf`）。
2. 运行流程：

```bash
bash pipe/1-run.sh conf/1-run.conf
```

如需强制重跑，将配置中的 `FORCE` 设为 `true`。

## 配置项说明（1-run.conf）

`1-run.sh` 通过 `source` 读取配置，常用字段如下：

- `PROJECT_DIR`：模块根目录。
- `VCF_PATH`：输入 VCF（可为 bgzip 压缩）。
- `SAMPLE_TABLE_PATH`：样本信息表。
- `FILTERED_SAMPLE_TABLE_PATH`：过滤后的样本表输出路径。
- `OUTPUT_DIR`：主输出目录。
- `LOG_DIR`：日志目录。
- `GROUP_COLS`：分群列名，逗号分隔（例如 `Population,Region`）。
- `ID_COL`：样本 ID 列名。
- `CHROM_NAME`：染色体名（如 `chrM`）。
- `USE_FILTERED_TABLE`：是否执行样本表过滤（`true/false`）。
- `SKIP_COMPLEX_VARIANTS`：是否跳过复杂变异（`true/false`）。
- `ENABLE_SIGNIFICANCE`：是否启用显著性分析（`true/false`）。
- `PERMUTATION_N`：置换次数（默认 1000）。
- `BOOTSTRAP_N`：bootstrap 次数（默认 1000）。
- `RANDOM_SEED`：随机种子（可选）。
- `CONDA_ENV`：运行脚本的 conda 环境名。
- `FORCE`：是否忽略成功日志强制重跑（`true/false`）。

## 执行流程说明

1. 初始化日志并记录参数。
2. 若日志中已出现 `SUCCESS` 且 `FORCE != true`，则跳过执行。
3. 可选：运行 `python/filter_sample_table.py` 过滤样本表。
4. 运行主分析脚本 `python/population_genetics.py`。
5. 完成后写入 `SUCCESS` 标记。

## 日志与断点续跑

- 日志文件：`log/1-run.sh.log`（自动追加）。
- 断点续跑：检测到 `SUCCESS` 将自动退出。

## 常见注意事项

- 每个群体样本数需 ≥2，否则主程序会报错。
- 变异过滤策略会显著影响统计量，建议固定过滤规则后再做比较。
```
