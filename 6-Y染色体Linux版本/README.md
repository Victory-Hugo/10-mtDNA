# Haplogroup Population Profile Generator

## 简介
单倍群群体频率分析工具 - Linux模块化版本  

## 快速使用

### 1. 准备数据
将您的数据文件放在 `input/` 目录中，支持的格式：
- `Table_Haplogroup.xlsx` (推荐)
- `Table_Haplogroup.xls`
- `Table_Haplogroup.csv`
- `Table_Haplogroup.txt`

**数据格式要求**：必须包含以下三列（列名不区分大小写）：
- `sample_id` / `sample` / `id` - 样本ID
- `haplogroup_full` / `haplogroup` / `hap` - 单倍群名称
- `population` / `pop` / `group` - 群体名称

### 2. 运行分析
```bash
./run_haplo_app.sh
```

### 3. 查看结果
结果将保存在 `output/` 目录中，包含：
- `output/out_NC01/` - NC=1的分析结果
- `output/out_NC02/` - NC=2的分析结果
- `output/out_NC03/` - NC=3的分析结果

每个目录包含：
- `tables/` - 数据表格（频率矩阵、距离矩阵等）
- `figs/` - 图表（热图、系统树、MDS图等）

## 系统要求
- Linux系统
- R >= 3.5.0
- 必要的R包（脚本会自动检查和安装）

## 文件说明
- `1-utils.R` - 通用工具函数
- `2-data_preprocessing.R` - 数据预处理
- `3-frequency_calculation.R` - 频率计算
- `4-distance_calculation.R` - 距离计算
- `5-visualization.R` - 数据可视化
- `run_haplo_app.sh` - 主控制脚本

## 配置
可编辑 `run_haplo_app.sh` 中的配置部分（约第175行）：
```bash
readonly NC_ARRAY=(1 2 3)  # NC范围
readonly MAX_PARALLEL=4    # 并行数
```

---
© 2025 BigLin. 基于单倍群的群体遗传学分析工具。