# -----------------------------------------------------------------------------
# 脚本功能说明:
# 本脚本用于调用指定的 Python 解释器和脚本，将 Excel 文件中的数据转换为 ID_Hap 格式，并输出到指定的文本文件。
#
# 参数说明:
#   PYTHON3    - Python 解释器路径（可根据实际环境修改）
#   PYTHON_SRC - Python 源代码路径（实现 Excel 到 ID_Hap 的转换）
#   EXCEL_PATH - 输入的 Excel 文件路径
#   OUT        - 输出的文本文件路径
#
# 使用方法:
#   直接运行本脚本即可完成数据转换，无需手动传递参数。
#
# 依赖:
#   需要安装 pandas 库（可通过 conda 安装）
#
# 注意事项:
#   请根据实际环境修改 Python 解释器和脚本路径。
# -----------------------------------------------------------------------------
#!/bin/bash
#TODO conda install pandas

PYTHON3="/home/luolintao/miniconda3/envs/pyg/bin/python3" #todo 替换为python解释器
PYTHON_SRC="/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/2-母系专用/python/python：Excel→ID_Hap.py"
EXCEL_PATH="/mnt/f/OneDrive/文档（共享）/2_mtDNA/线粒体_现代DNA83.xlsx"
OUT="/mnt/f/4_20K_CPGDP/0-DAPC/conf/ID_Hap.txt"

"$PYTHON3" \
    "$PYTHON_SRC" \
    "$EXCEL_PATH" \
    "$OUT"

