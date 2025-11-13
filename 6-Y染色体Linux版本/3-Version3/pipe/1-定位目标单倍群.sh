: '
脚本名称: 1-定位目标单倍群.sh

功能简介:
本脚本用于调用指定的 Python 脚本，对线粒体单倍群数据进行整理和定位，输出目标单倍群相关的结果文件。主要处理流程包括单倍群名称订正、等级排序、未查询到和无法处理的样本统计等。

使用说明:
1. 请根据实际环境修改 PYTHON3、SRC_DIR、ID_HAP_TXT、OUT_DIR 等变量为正确的路径。
2. 输入文件包括样本单倍群信息、订正表、phylotree、目标单倍群列表和错误纠正表。
3. 输出文件包括订正后的单倍群名称、正序/逆序等级、未查询到的样本、最终结果和无法处理的样本。
4. 脚本分两次运行，分别针对不同的目标单倍群配置（Type1 和 Type2）。

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
# todo conda install pandas openpyxl

# ========== Python 配置 ==========
PYTHON3="/home/luolintao/miniconda3/envs/BigLin/bin/python3" #todo 替换为python解释器
SRC_DIR="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/6-Y染色体Linux版本/3-Version3/"  #todo 替换为/2-母系专用 绝对路径
EXCEL_FILE="${SRC_DIR}/input/Table_Haplogroup.xlsx" #todo Excel文件路径,第一列是ID，第二列是Haplogroup，第三列是Population
OUT_DIR="${SRC_DIR}/output/" #todo 输出文件路径

# ========== R 配置 ==========
RSCRIPT="/home/luolintao/miniconda3/envs/BigLin/bin/Rscript" #todo 替换为Rscript路径
R_SCRIPT="${SRC_DIR}/pipe/1-主执行函数_MTDNA.R" #todo R脚本路径（参数化版本）
META_CSV="${SRC_DIR}/conf/Group.csv" #todo 群体元数据文件（可选）
YSTR_CSV="${SRC_DIR}/input/YSTR_table.csv" #todo Y-STR数据文件（可选）


#* ========== 自动生成 ============
ID_HAP_TXT="${SRC_DIR}/input/ID_Hap.txt"
HAPLOGROUP_HIERARCHY_CSV="${SRC_DIR}/input/Example_Haplogroup_Hierarchy.csv"


# *从Excel提取ID和Haplogroup列
echo "正在从Excel文件提取数据..."
"${PYTHON3}" "${SRC_DIR}/python/0-extract_from_excel.py" \
  --input-excel "${EXCEL_FILE}" \
  --output-txt "${ID_HAP_TXT}"

if [ $? -ne 0 ]; then
    echo "错误: 从Excel提取数据失败"
    exit 1
fi

# *请勿修改下列内容
"${PYTHON3}" "${SRC_DIR}/python/1-Seek_levels.py" \
  --file-id-hap "${ID_HAP_TXT}" \
  --file-correction "${SRC_DIR}/conf/Haplogrep单倍群分型订正表.csv" \
  --file-phylotree "${SRC_DIR}/conf/线粒体单倍群phylotree(version17)2025年3月12日.txt" \
  --file-target "${SRC_DIR}/conf/目标_Type1.txt" \
  --file-rough-fix "${SRC_DIR}/conf/错误纠正_Type1.txt" \
  --out-fixed-hap "${OUT_DIR}/订正之后的单倍群名称.txt" \
  --out-forward-level "${OUT_DIR}/正序等级.txt" \
  --out-reversed-level "${OUT_DIR}/逆序等级.txt" \
  --out-not-found "${OUT_DIR}/没有查询到请核实.txt" \
  --out-final "${OUT_DIR}/最终_Type1.txt" \
  --out-unmatched "${OUT_DIR}/无法处理.txt"

# *请勿修改下列内容
"${PYTHON3}" "${SRC_DIR}/python/1-Seek_levels.py" \
  --file-id-hap "${ID_HAP_TXT}" \
  --file-correction "${SRC_DIR}/conf/Haplogrep单倍群分型订正表.csv" \
  --file-phylotree "${SRC_DIR}/conf/线粒体单倍群phylotree(version17)2025年3月12日.txt" \
  --file-target "${SRC_DIR}/conf/目标_Type2.txt" \
  --file-rough-fix "${SRC_DIR}/conf/错误纠正_Type2.txt" \
  --out-fixed-hap "${OUT_DIR}/订正之后的单倍群名称.txt" \
  --out-forward-level "${OUT_DIR}/正序等级.txt" \
  --out-reversed-level "${OUT_DIR}/逆序等级.txt" \
  --out-not-found "${OUT_DIR}/没有查询到请核实.txt" \
  --out-final "${OUT_DIR}/最终_Type2.txt" \
  --out-unmatched "${OUT_DIR}/无法处理.txt"

# *补齐至下游单倍群
echo "正在补齐至下游单倍群..."
"${PYTHON3}" "${SRC_DIR}/python/2-fill_down.py" \
  --input-file "${OUT_DIR}/正序等级.txt" \
  --output-file "${OUT_DIR}/正序等级-补齐至下游.txt"

# *按首字母调整等级
echo "正在按首字母调整等级..."
"${PYTHON3}" "${SRC_DIR}/python/3-adjust_header.py" \
  --id-hap-file "${ID_HAP_TXT}" \
  --forward-level-file "${OUT_DIR}/正序等级-补齐至下游.txt" \
  --output-file "${HAPLOGROUP_HIERARCHY_CSV}"

echo "删除不必要的文件……"
rm -rf "${OUT_DIR}"
echo "Python流程运行完毕 ✅"

# ========== 运行R分析 ==========
echo "开始运行R分析..."

# 调用R脚本并传递参数
"${RSCRIPT}" "${R_SCRIPT}" \
  "wd=${SRC_DIR}/pipe" \
  "input_xlsx=${EXCEL_FILE}" \
  "meta_csv=${META_CSV}" \
  "ystr_csv=${YSTR_CSV}" \
  "haplogroup_hierarchy_csv=${HAPLOGROUP_HIERARCHY_CSV}" \
  "root_out=${SRC_DIR}/output"

