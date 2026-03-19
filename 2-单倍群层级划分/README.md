# mtDNA 单倍群层级划分 Pipeline

该项目用于对人类 mtDNA 单倍群结果进行标准化、映射到 PhyloTree Build 17，并进一步输出 LLT 与 YuChunLi 两套宏单倍群结果。

## 功能

- 读取输入的 `tsv` 或 `csv` 文件
- 自动识别是否存在表头
- 使用 `data/Haplogrep单倍群分型订正表.tsv` 纠正输入单倍群名称
- 按 `data/phylotree_build17.tsv` 进行层级挂接
- 输出：
  - 原始单倍群
  - 订正后单倍群
  - PhyloTree 单倍群
  - LLT 宏单倍群
  - YuChunLi 宏单倍群
  - PhyloTree 级别、父节点、定义变异

## 目录结构

```text
2-单倍群层级划分/
├── conf/
│   └── Config.yaml
├── data/
│   ├── phylotree_build17.tsv
│   ├── phylotree_index.json
│   ├── Haplogrep单倍群分型订正表.tsv
│   ├── 目标_LLT.tsv
│   ├── 错误纠正_LLT.tsv
│   ├── 目标_YuChunLi.tsv
│   └── 错误纠正_YuChunLi.tsv
├── input/
├── output/
├── pipe/
│   └── run_pipeline.sh
├── python/
│   ├── build_tree_index.py
│   └── annotate_haplogroups.py
├── script/
│   ├── check_env.sh
│   └── load_config.sh
└── src/
```

## 主要模块

- `python/build_tree_index.py`
  - 双模式模块
  - 用于构建可复用的 `data/phylotree_index.json`
- `python/annotate_haplogroups.py`
  - 双模式模块
  - 用于读取样本输入并输出注释结果
- `pipe/run_pipeline.sh`
  - 总控脚本
  - 从 `conf/Config.yaml` 读取配置并调用两个 Python 模块

## 输入格式

输入文件默认是 `input/ID_Hap.tsv`，要求前两列为：

1. 样本 ID
2. 单倍群名称

支持：

- `tsv`
- `csv`
- 有表头
- 无表头

## 输出文件

默认输出文件：

- `output/mtDNA_haplogroup_annotation.tsv`

默认索引文件：

- `data/phylotree_index.json`

输出列包括：

- `ID`
- `Haplogroup_Original`
- `Haplogroup_Standardized`
- `Haplogroup_PhyloTree`
- `Haplogroup_LLT`
- `Haplogroup_YuChunLi`
- `Resolution_Status`
- `Phylotree_Level`
- `Phylotree_Parent`
- `Phylotree_Mutations`

## 运行方式

在项目根目录执行：

```bash
bash pipe/run_pipeline.sh
```

也可以显式指定配置文件：

```bash
bash pipe/run_pipeline.sh --config conf/Config.yaml
```

## 配置说明

配置文件位于 `conf/Config.yaml`。当前默认项包括：

- `paths.input_file`
- `paths.output_dir`
- `paths.tree_index_file`
- `paths.phylotree_table`
- `paths.input_correction`
- `paths.llt_targets`
- `paths.llt_correction`
- `paths.yuchunli_targets`
- `paths.yuchunli_correction`
- `runtime.output_prefix`
- `runtime.rebuild_index`
- `runtime.log_level`

其中：

- 当 `runtime.rebuild_index: true` 时，每次运行都会重建 `data/phylotree_index.json`
- 当 `runtime.rebuild_index: false` 时，会优先复用已有索引

## 单独运行模块

构建索引：

```bash
python3 python/build_tree_index.py \
  --phylotree-table data/phylotree_build17.tsv \
  --input-correction data/Haplogrep单倍群分型订正表.tsv \
  --llt-targets data/目标_LLT.tsv \
  --llt-correction data/错误纠正_LLT.tsv \
  --yuchunli-targets data/目标_YuChunLi.tsv \
  --yuchunli-correction data/错误纠正_YuChunLi.tsv \
  --output data/phylotree_index.json
```

样本注释：

```bash
python3 python/annotate_haplogroups.py \
  --input input/ID_Hap.tsv \
  --tree-index data/phylotree_index.json \
  --output output/mtDNA_haplogroup_annotation.tsv
```

## 备注

- 旧版脚本仍保留在 `python/` 中，当前流程未直接依赖它们。
- `phylotree_index.json` 放在 `data/` 下，方便后续复用和其他脚本直接调用。
