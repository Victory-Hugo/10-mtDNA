# 1-run-pipe.sh 使用说明

## 脚本简介

本脚本为“样本筛选与单倍群分析流程”的主控管道脚本，自动化完成从原始数据到频率分析、PCA降维及可视化的全流程。用户只需配置参数并运行本脚本，即可依次完成各分析步骤。

## 使用方法

**编辑参数**：所有用户参数已外部化，存放于 `conf/config.cmf` 文件中。请根据实际需求编辑该配置文件。

## 参数说明（均在 conf/config.cmf 中配置）

- `BASE_DIR`：包含完整的基本信息表。
- `PYTHON3`：指定 Python3 解释器路径。
- `INPUT_EXCEL`：输入的 Excel 文件（可为相对路径或绝对路径）。
- `FORCE_RUN`：是否强制重新运行所有步骤（true/false）。
- `NEW_SAMPLE_SOURCE`：新样本来源（用于样本筛选）。
- `COUNTRY_EXTRA_COLS`、`ETHNICITY_EXTRA_COLS`：国家/民族附加列名（数组）。
- `SELECTED_COLUMNS`：频率分析所需的列名（数组）。
- `MIN_SAMPLE_COUNT`：频率分析的最小样本数阈值。
- `HAPLOGROUP_COLUMN`：单倍群列名。
- `HEATMAP_GROUP_COL`：热图分组列名。
- `CLASS_BIG_COLUMN`：PCA后处理的大类列名。
- `VISUAL_COLOR_CSV`：频率柱状图颜色映射文件。
- `VISUAL_GROUP_COLUMN`：频率柱状图分组列名。
- `PCA_VIZ_COLOR_CSV`：PCA可视化颜色映射文件。
- `RSCRIPT`：Rscript 命令路径。
- 其他参数详见 `conf/config.cmf` 注释。

## 脚本功能与流程

脚本共分为 12 个主要步骤，自动串联执行：

1. **单倍群整理**：调用 `haplogroup_organizer.py`，整理单倍群信息。
   先标准化单倍群的分型结果，再合并LLT和LYC两个分级下的宏单倍群。
2. **样本筛选和清洗**：调用 `data_cleaner.py`，筛选并清洗样本。
   （1）提取新样本以及不是由LLT开头的公共样本
   （2）标准化这些样本的基础信息，比如统一国家、省份、民族的表达方式
   （3）从订正表中提取需要的列（南北、大洲、语系等）
   （4）将未知的信息改为`Unknown`
   （5）Classification:中国的用population_province分组，外国的用population_country分组
3. **合并单倍群信息**：调用 `merger.py`，合并样本与单倍群信息。
4. **频率分析前置处理**：调用 `frequency_preprocessor.py`，生成频率矩阵及分组映射。
   （1）删除有未知信息的样本（对于基于频率的分析没有帮助）
   （2）创建一个Group列，其中汉语族为Sinitic_North或者Sinitic_South，其他语言人群则以他们的语系为分组（常用的分组方式）。
   （3）删除人群（Classification）数小于20的人群（可修改，建议≥20）
   （4）以人群（Classification）为单位，计算LYC（可修改）单倍群的频率
5. **PCA降维（全球）**：调用 `pca_core.py`，对全球数据做PCA降维。
6. **PCA后处理（全球）**：调用 `pca_postprocessor.py`，整理PCA结果。
7. **PCA降维（中国）**：对中国数据做PCA降维。
8. **PCA结果处理（中国）**：整理中国PCA结果。
   整理降维结果，为后续的PCA可视化做准备：Class_small为可视化设置不同图形的依据，这里Classification为Class_small；Class_big为可视化时着色的依据，这里以Group作为着色依据（可修改）
9. **频率堆叠柱状图（全球+中国）**：调用 `haplogroup_barplot.py`，生成频率堆叠柱状图。需要输入以Haplogroup和Color为列名的颜色csv文件（必须包含这两列）。
10. **PCA可视化（R脚本）**：调用 `pca_visualization.R`，生成PCA可视化图。这里需要一个以Class_big,color为列名的颜色csv文件（必须包含这两列）。
11. **频率热图可视化（R脚本）**：调用 `heatmap_visualization.R`，生成频率热图。
12. **四川省各市样本数量柱状图**：调用 `sichuan_city_barplot.py`，以第2步输出的样本基本信息为输入，筛选四川省样本，按城市统计样本数量，生成横轴为 City、纵轴为 Sichuan_Count 的柱状图。输出 PDF 图和统计 CSV，保存至 `output/四川省各市样本数量/`。

每一步均有日志记录，支持断点续跑和强制重跑。

## 日志与错误处理

- 日志文件：`log/run-pipe.sh.log`
- 每步执行情况、错误信息均会写入日志。
- 若遇到缺失文件或执行失败，脚本会自动终止并输出错误信息。

## 依赖环境

- Python 3 及相关依赖包
- R 及 Rscript 命令、相关R包

## 常见问题

- 配置文件缺失或参数错误会导致脚本无法启动。
- 各分析脚本需具备可执行权限，且依赖包需提前安装。
- 如需跳过已完成步骤，设置 `FORCE_RUN=false`；如需全流程重跑，设为 `true`。

---
如需详细参数说明或遇到问题，请查阅 `conf/config.cmf` 注释或联系维护者。