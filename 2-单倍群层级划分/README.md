# 🧬 mtDNA 单倍群层级划分 Pipeline

![Landscape.png](https://picturerealm.oss-cn-chengdu.aliyuncs.com/obsidian/20260319171335158.png)

一个面向**人类 mtDNA 单倍群注释与层级归类**的生物信息学流程。  
该工具可将输入的单倍群名称进行标准化，映射到 **PhyloTree Build 17**，并进一步输出两套宏单倍群结果：

- 🏷️ **LLT**
- 🏷️ **YuChunLi**

---

## ✨ 功能概览

本流程支持以下任务：

- 📥 读取输入文件，支持 `tsv` 与 `csv`
- 🧾 自动识别是否存在表头
- 🛠️ 使用订正表纠正不规范单倍群名称
- 🌳 将单倍群挂接到 **PhyloTree Build 17**
- 🧬 输出标准化后的细粒度单倍群信息
- 🧭 同时映射为 **LLT** 与 **YuChunLi** 两套宏单倍群
- 📊 输出层级、父节点、定义变异等注释信息

---

## 🧱 目录结构

```text
2-单倍群层级划分
├── conf
│   └── Config.yaml
├── data
│   ├── phylotree_build17.tsv
│   ├── phylotree_index.json
│   ├── Haplogrep单倍群分型订正表.tsv
│   ├── 目标_LLT.tsv
│   ├── 错误纠正_LLT.tsv
│   ├── 目标_YuChunLi.tsv
│   └── 错误纠正_YuChunLi.tsv
├── input
│   └── ID_Hap.tsv
├── output
│   └── mtDNA_haplogroup_annotation.tsv
├── pipe
│   └── run_pipeline.sh
├── python
│   ├── build_tree_index.py
│   ├── annotate_haplogroups.py
│   ├── merge_new_haplogroups.py
│   └── parse_phylotree.py
├── script
│   ├── check_env.sh
│   ├── console_ui.sh
│   └── load_config.sh
└── README.md
````

---

## 🚀 流程逻辑

整个流程分为两个核心 Python 模块：

### 1️⃣ `build_tree_index.py`

用于构建可复用的树索引文件 `phylotree_index.json`。

主要完成：

* 读取 PhyloTree 表
* 读取输入单倍群订正表
* 读取 LLT 目标集合与纠错表
* 读取 YuChunLi 目标集合与纠错表
* 构建可快速查询的层级索引

### 2️⃣ `annotate_haplogroups.py`

用于读取样本输入文件，并输出单倍群注释结果。

主要完成：

* 读取输入样本
* 自动识别表头
* 标准化输入单倍群
* 查询树索引
* 输出 PhyloTree 注释结果
* 输出 LLT 与 YuChunLi 宏单倍群结果

---

## 📥 输入说明

默认输入文件为：

```text
input/ID_Hap.tsv
```

要求前两列分别为：

1. 样本 ID
2. 单倍群名称

支持以下格式：

* `tsv`
* `csv`
* 有表头
* 无表头

---

## 🛠️ 数据文件说明

### 🌳 `phylotree_build17.tsv`

PhyloTree 主表。
包含单倍群名称、层级、父节点、定义变异等信息，是整个层级系统的核心参考。

### 🩹 `Haplogrep单倍群分型订正表.tsv`

输入纠正表。
第一列为可能出现的不标准名称，第二列为标准单倍群名称。
用于将输入单倍群统一映射为标准命名。

### 🎯 `目标_LLT.tsv`

LLT 套餐中的目标宏单倍群列表。
每一行是一个允许输出的宏单倍群名称。

### 🔁 `错误纠正_LLT.tsv`

LLT 纠错表。
第一列为标准单倍群名称，第二列为需要统一映射到的 LLT 宏单倍群名称。

### 🎯 `目标_YuChunLi.tsv`

YuChunLi 套餐中的目标宏单倍群列表。
每一行是一个允许输出的宏单倍群名称。

### 🔁 `错误纠正_YuChunLi.tsv`

YuChunLi 纠错表。
第一列为标准单倍群名称，第二列为需要统一映射到的 YuChunLi 宏单倍群名称。

---

## 🔍 注释逻辑说明

对每个输入单倍群，流程按以下顺序处理：

### 第一步：输入名称标准化

先用 `Haplogrep单倍群分型订正表.tsv` 对输入名称进行纠正，得到标准单倍群名称。

### 第二步：挂接到 PhyloTree

将标准化后的单倍群映射到 `phylotree_build17.tsv` 中，获取：

* 单倍群名称
* level
* parent
* mutations

### 第三步：映射到 LLT

沿着该单倍群的祖先路径向上回溯，寻找第一个属于 `目标_LLT.tsv` 的节点。
若该节点命中 `错误纠正_LLT.tsv`，则按纠错表替换为最终 LLT 宏单倍群。

### 第四步：映射到 YuChunLi

沿着该单倍群的祖先路径向上回溯，寻找第一个属于 `目标_YuChunLi.tsv` 的节点。
若该节点命中 `错误纠正_YuChunLi.tsv`，则按纠错表替换为最终 YuChunLi 宏单倍群。

---

## 📤 输出说明

默认输出文件为：

```text
output/mtDNA_haplogroup_annotation.tsv
```

输出列包括：

* `ID`
* `Haplogroup_Original`
* `Haplogroup_Standardized`
* `Haplogroup_PhyloTree`
* `Haplogroup_LLT`
* `Haplogroup_YuChunLi`
* `Resolution_Status`
* `Phylotree_Level`
* `Phylotree_Parent`
* `Phylotree_Mutations`

---

## 📌 `Resolution_Status` 含义

该字段用于标识输入单倍群的解析状态，便于后续质控与排查。

常见状态包括：

* ✅ `exact`
  输入或订正后的单倍群可直接在树中找到

* 🩹 `corrected`
  输入单倍群先经过订正，再成功映射到树中

* 🌿 `inferred_ancestor`
  未能直接匹配时，向上回退到可识别祖先节点后完成注释

* ❓ `unresolved`
  无法完成有效映射

---

## ⚙️ 配置文件

配置文件位于：

```text
conf/Config.yaml
```

当前主要配置项包括：

### `tools`

* `python`
  Python 解释器名称

### `paths`

* `input_file`
* `output_dir`
* `tree_index_file`
* `phylotree_table`
* `input_correction`
* `llt_targets`
* `llt_correction`
* `yuchunli_targets`
* `yuchunli_correction`

### `runtime`

* `output_prefix`
* `rebuild_index`
* `log_level`

---

## ▶️ 运行方式

在项目根目录执行：

```bash
bash pipe/run_pipeline.sh
```

也可以显式指定配置文件：

```bash
bash pipe/run_pipeline.sh --config conf/Config.yaml
```

---

## 🧪 单独运行模块

### 构建树索引

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

### 样本单倍群注释

```bash
python3 python/annotate_haplogroups.py \
  --input input/ID_Hap.tsv \
  --tree-index data/phylotree_index.json \
  --output output/mtDNA_haplogroup_annotation.tsv
```

---

## 💡 设计特点

* 🧩 **双模块设计**
  建索引与样本注释分离，便于维护与扩展

* ♻️ **索引可复用**
  树索引可重复使用，避免每次都重新解析全表

* 🔧 **配置驱动**
  关键路径写入配置文件，便于迁移与部署

* 📚 **适合继续扩展**
  后续可加入更多宏单倍群体系、统计报告、质控文件与日志输出

---

## 📝 备注

* 当前流程以 **人类 mtDNA 单倍群** 为目标对象
* `phylotree_index.json` 放在 `data` 中，便于复用
* `python` 目录中保留部分早期脚本，当前主流程依赖的是：

  * `build_tree_index.py`
  * `annotate_haplogroups.py`

---

## ✅ 总结

这个 Pipeline 的核心目标是：

**将输入样本的 mtDNA 单倍群名称标准化，并稳定映射到 PhyloTree 与两套宏单倍群体系，输出结构化注释结果。**

适合用于：

* 人群遗传学分析
* 母系谱系分类整理
* 宏单倍群统一输出
* 下游统计分析前的数据标准化

🧬📦🌳

