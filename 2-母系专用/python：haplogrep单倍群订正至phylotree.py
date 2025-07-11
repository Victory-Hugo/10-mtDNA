import pandas as pd

# 加载订正表
correction_table_path = 'F:/OneDrive/文档（科研）/脚本/我的科研脚本/Python/母系专用/Haplogrep单倍群分型订正表.csv'
correction_table = pd.read_csv(correction_table_path)

# 加载主要数据文件,数据中需要修改的列名是 'Haplogroup'
input_file_path = 'C:/Users/victo/Desktop/新建 Text Document.txt'
data = pd.read_csv(input_file_path, sep='\t',encoding="UTF-8")

# 创建字典以进行查找和替换
correction_dict = dict(zip(correction_table['Origin'], correction_table['Revised']))

# 替换主要数据文件中 'Haplogroup' 列的值
data['Haplogroup'] = data['Haplogroup'].replace(correction_dict)
# 打印被替换的具体值
print(data[data['Haplogroup'].isin(correction_dict.values())])

# 保存修改后的数据文件
output_file_path = input_file_path
data.to_csv(output_file_path, sep='\t', index=False)
