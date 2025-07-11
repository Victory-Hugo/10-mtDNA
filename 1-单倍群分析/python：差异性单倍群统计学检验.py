import pandas as pd
from scipy.stats import chi2_contingency, fisher_exact

# 读取数据
data = pd.read_csv('C:/Users/victo/Desktop/Haplogroup_counts.txt', sep='\t')

# 从Category_totals.txt中读取总和数据
category_totals = pd.read_csv('C:/Users/victo/Desktop/Category_totals.txt', sep='\t', header=None, index_col=0)
category_totals.columns = ['Total']

# 提供列名供用户选择
columns = data.columns[1:]  # 排除第一列（单倍群名称）

# 显示可选列名
print("可用的列名如下：")
for i, col in enumerate(columns, start=1):
    print(f"{i}: {col}")

# 交互式选择列
col1_index = int(input("请选择第一列的编号：")) - 1
col2_index = int(input("请选择第二列的编号：")) - 1

col1_name = columns[col1_index]
col2_name = columns[col2_index]

# 获取外部总和数据
total_col1 = category_totals.loc[col1_name, 'Total']
total_col2 = category_totals.loc[col2_name, 'Total']

# 计算每个单倍群在所选群体中的频率
data[f'{col1_name}_freq'] = data[col1_name] / total_col1
data[f'{col2_name}_freq'] = data[col2_name] / total_col2

# 计算频率差异
data['freq_diff'] = abs(data[f'{col1_name}_freq'] - data[f'{col2_name}_freq'])

# 找出频率差异最大的前十个单倍群
rank_number = int(input("请输入需要计算的前几名差异单倍群：（例如输入10）"))
top_10_diff = data.nlargest(rank_number, 'freq_diff')

# 进行卡方检验
top_10_diff['chi2_p_value'] = top_10_diff.apply(
    lambda row: chi2_contingency([[row[col1_name], total_col1 - row[col1_name]], 
                                  [row[col2_name], total_col2 - row[col2_name]]])[1],
    axis=1
)

# 进行Fisher精确检验
top_10_diff['fisher_p_value'] = top_10_diff.apply(
    lambda row: fisher_exact([[row[col1_name], total_col1 - row[col1_name]], 
                              [row[col2_name], total_col2 - row[col2_name]]])[1],
    axis=1
)

# 保存结果至桌面的TXT文件
output_path = 'C:/Users/victo/Desktop/haplogroup_difference_result.txt'
top_10_diff.to_csv(output_path, sep='\t', index=False)

print(f"结果已保存至 {output_path}")
