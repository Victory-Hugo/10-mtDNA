: '
脚本名称: 1-定位目标单倍群.sh

功能简介:
本脚本用于调用指定的 Python 脚本，对线粒体单倍群数据进行整理和定位，输出目标单倍群相关的结果文件。主要处理流程包括单倍群名称订正、等级排序、未查询到和无法处理的样本统计等。

使用说明:
1. 请根据实际环境修改 PYTHON3、SRC_DIR、ID_HAP_TXT、OUT_DIR 等变量为正确的路径。
2. 输入文件包括样本单倍群信息、订正表、phylotree、目标单倍群列表和错误纠正表。
3. 输出文件包括订正后的单倍群名称、正序/逆序等级、未查询到的样本、最终结果和无法处理的样本。
4. 脚本分两次运行，分别针对不同的目标单倍群配置（YuChunLi 和 LLT）。

参数说明:
- --file-id-hap: 样本单倍群信息文件，第一列为ID，第二列为单倍群全称，制表符分割。
- --file-correction: 单倍群分型订正表。
- --file-phylotree: 线粒体单倍群phylotree信息。
- --file-target: 目标单倍群列表。
- --file-rough-fix: 错误纠正表。
- --out-fixed-hap: 输出订正后的单倍群名称。
- --out-forward-level: 输出正序等级。
- --out-reversed-level: 输出逆序等级。
- --out-not-found: 输出未查询到的样本。
- --out-final: 输出最终结果。
- --out-unmatched: 输出无法处理的样本。

注意事项:
- 请勿修改脚本中标注为“请勿修改下列内容”的部分，以保证脚本正常运行。
- 需提前安装 pandas 库（可通过 conda install pandas）。
'
#!/bin/bash
# todo conda install pandas

PYTHON3="/home/luolintao/miniconda3/envs/pyg/bin/python3" #todo 替换为python解释器
SRC_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/2-母系专用/"  #todo 替换为/2-母系专用 绝对路径
ID_HAP_TXT="/mnt/f/6_起源地混合地/3-发育树可视化/data/ID_Hap.txt" # todo 第一列为ID，第二列为单倍群全称，列名为ID Haplogroup，制表符分割
OUT_DIR="/mnt/f/6_起源地混合地/3-发育树可视化/data/" #todo 输出文件路径


# *请勿修改下列内容
"${PYTHON3}" "${SRC_DIR}/python/python：快速整理线粒体单倍群至目标单倍群.py" \
  --file-id-hap "${ID_HAP_TXT}" \
  --file-correction "${SRC_DIR}/conf/Haplogrep单倍群分型订正表.csv" \
  --file-phylotree "${SRC_DIR}/conf/线粒体单倍群phylotree(version17)2025年3月12日.txt" \
  --file-target "${SRC_DIR}/conf/目标_YuChunLi.txt" \
  --file-rough-fix "${SRC_DIR}/conf/错误纠正_YuChunLi.txt" \
  --out-fixed-hap "${OUT_DIR}/订正之后的单倍群名称.txt" \
  --out-forward-level "${OUT_DIR}/正序等级.txt" \
  --out-reversed-level "${OUT_DIR}/逆序等级.txt" \
  --out-not-found "${OUT_DIR}/没有查询到请核实.txt" \
  --out-final "${OUT_DIR}/最终_YuChunLi.txt" \
  --out-unmatched "${OUT_DIR}/无法处理.txt"

# *请勿修改下列内容
"${PYTHON3}" "${SRC_DIR}/python/python：快速整理线粒体单倍群至目标单倍群.py" \
  --file-id-hap "${ID_HAP_TXT}" \
  --file-correction "${SRC_DIR}/conf/Haplogrep单倍群分型订正表.csv" \
  --file-phylotree "${SRC_DIR}/conf/线粒体单倍群phylotree(version17)2025年3月12日.txt" \
  --file-target "${SRC_DIR}/conf/目标_LLT.txt" \
  --file-rough-fix "${SRC_DIR}/conf/错误纠正_LLT.txt" \
  --out-fixed-hap "${OUT_DIR}/订正之后的单倍群名称.txt" \
  --out-forward-level "${OUT_DIR}/正序等级.txt" \
  --out-reversed-level "${OUT_DIR}/逆序等级.txt" \
  --out-not-found "${OUT_DIR}/没有查询到请核实.txt" \
  --out-final "${OUT_DIR}/最终_LLT.txt" \
  --out-unmatched "${OUT_DIR}/无法处理.txt"

