import pandas as pd

# 加载数据
file_path = 'C:/Users/victo/Desktop/线粒体_现代DNA.csv'  # 请替换为你的文件路径
data = pd.read_csv(file_path)

# 定义单倍群列名
haplogroup_column = "HAPLOGROUP(17.2)"

# 按单倍群进行50%的抽样
sampled_data = data.groupby(haplogroup_column).apply(lambda x: x.sample(frac=0.5, random_state=1))
sampled_data = sampled_data.droplevel(0)  # 移除多余的索引

# 过滤掉以 "L0", "L1", "L2" 开头的单倍群
filtered_sampled_data = sampled_data[~sampled_data[haplogroup_column].str.startswith(('L0', 'L1', 'L2'))]

# 保存结果为CSV文件
output_path = 'C:/Users/victo/Desktop/filtered_sampled_data.csv'  # 请替换为你想保存的路径
filtered_sampled_data.to_csv(output_path, index=False)

# 输出剩余样本数量
print(f"过滤后的样本数量: {filtered_sampled_data.shape[0]}")
