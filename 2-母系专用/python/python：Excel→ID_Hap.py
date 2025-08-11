import pandas as pd
import sys


# 从命令行读取参数
if len(sys.argv) != 3:
    print(f"用法: python {sys.argv[0]} <输入Excel路径> <输出文件路径>")
    sys.exit(1)

input_excel = sys.argv[1]  # 输入 Excel 文件
output_file = sys.argv[2]  # 输出文件路径

# 读取数据
df_META = pd.read_excel(input_excel, sheet_name="现代线粒体DNA汇总")
df_QC = pd.read_excel(input_excel, sheet_name="质量控制")
df_Chip = pd.read_excel(input_excel, sheet_name="芯片单倍群")

df_META = df_META.loc[:, ['ID', 'Sequence_method']]
df_QC = df_QC.loc[:, ['SampleID', 'Haplogroup', 'Quality']]
df_Chip = df_Chip.loc[:, ['SampleID', 'Haplogroup', 'Quality']]

df_merge = pd.concat([df_Chip, df_QC], axis=0).dropna()
df_remerge = pd.merge(df_merge, df_META, left_on='SampleID', right_on='ID', how='outer')
df_remerge = df_remerge.dropna(subset=['ID'])

df_remerge.loc[:, ['ID', 'Haplogroup']].to_csv(output_file, sep='\t', index=False)
