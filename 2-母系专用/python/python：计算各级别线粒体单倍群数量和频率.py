import pandas as pd
import numpy as np
from multiprocessing import Pool, cpu_count

# 加载数据集
file_path = 'C:/Users/victo/Desktop/新建 Text Document.txt'  # 替换为你的数据文件路径
data = pd.read_csv(file_path,sep="\t")

# 列出单倍群字段
haplogroup_columns = [f'Level_{i}' for i in range(20)]  # Level_0到Level_19

# 获取所有单倍群
def get_unique_haplogroups(df, columns):
    unique_haplogroups = set()
    for col in columns:
        unique_haplogroups.update(df[col].dropna().unique())  # 将各列中的非空唯一值加入集合
    return unique_haplogroups

# 单个省份的统计工作（并行执行的任务）
def calculate_haplogroup_for_province(args):
    province, province_data, haplogroups = args
    province_haplogroup_data = province_data[haplogroup_columns].to_numpy()

    haplogroup_counts = {}
    for haplogroup in haplogroups:
        # 使用 NumPy 的向量化操作来检查单倍群是否出现在该省份的某行中
        mask = np.any(province_haplogroup_data == haplogroup, axis=1)
        # 统计出现次数（每一行最多计数一次）
        haplogroup_counts[haplogroup] = np.sum(mask)
    
    total_samples = len(province_data)
    haplogroup_frequencies = {haplogroup: count / total_samples for haplogroup, count in haplogroup_counts.items()}
    
    return province, haplogroup_counts, haplogroup_frequencies

# 多核计算函数
def calculate_haplogroup_counts_parallel(df, province_col):
    # 获取所有单倍群
    haplogroups = get_unique_haplogroups(df, haplogroup_columns)
    
    # 按省份分组数据
    province_groups = [(province, df[df[province_col] == province], haplogroups) for province in df[province_col].unique()]

    # 使用多核CPU进行并行计算
    with Pool(cpu_count()) as pool:
        results = pool.map(calculate_haplogroup_for_province, province_groups)

    # 汇总结果
    haplogroup_counts = {}
    haplogroup_frequencies = {}

    for province, counts, frequencies in results:
        haplogroup_counts[province] = counts
        haplogroup_frequencies[province] = frequencies

    return haplogroup_counts, haplogroup_frequencies

if __name__ == '__main__':
    # 交互式输入省份字段
    province_field = input("请输入要计算的字段名称，例如 Province_cn: ")

    # 计算结果
    haplogroup_counts, haplogroup_frequencies = calculate_haplogroup_counts_parallel(data, province_field)

    # 将数量结果转换为 DataFrame 以表格形式展示
    haplogroup_counts_df = pd.DataFrame(haplogroup_counts).T.fillna(0)

    # 将频率结果转换为 DataFrame 以表格形式展示
    haplogroup_frequencies_df = pd.DataFrame(haplogroup_frequencies).T.fillna(0)

    # 重置索引，使用户输入的字段成为第一列
    haplogroup_counts_df.index.name = province_field
    haplogroup_frequencies_df.index.name = province_field
    haplogroup_counts_df.reset_index(inplace=True)
    haplogroup_frequencies_df.reset_index(inplace=True)

    # 保存数量和频率表格为CSV文件
    haplogroup_counts_df.to_csv("C:/Users/victo/Desktop/haplogroup_counts.csv", index=False)
    haplogroup_frequencies_df.to_csv("C:/Users/victo/Desktop/haplogroup_frequencies.csv", index=False)

    # 输出数量和频率表格
    print("单倍群数量统计表格:")
    print(haplogroup_counts_df)
    print("\n单倍群频率表格:")
    print(haplogroup_frequencies_df)

#####以下是在普通计算机中使用的脚本################
# import pandas as pd
# import numpy as np

# # 加载数据集
# file_path = '你的数据路径.csv'  # 替换为你的数据文件路径
# data = pd.read_csv(file_path)

# # 列出单倍群字段
# haplogroup_columns = [f'Level_{i}' for i in range(20)]  # Level_0到Level_19

# # 获取所有单倍群
# def get_unique_haplogroups(df, columns):
#     unique_haplogroups = set()
#     for col in columns:
#         unique_haplogroups.update(df[col].dropna().unique())  # 将各列中的非空唯一值加入集合
#     return unique_haplogroups

# # 高效计算每个单倍群在每个省份的出现次数
# def calculate_haplogroup_counts(df, province_col):
#     # 创建一个字典来存储结果
#     haplogroup_counts = {}
    
#     # 获取所有单倍群
#     haplogroups = get_unique_haplogroups(df, haplogroup_columns)
    
#     # 对每个省份进行处理
#     for province in df[province_col].unique():
#         haplogroup_counts[province] = {}
#         # 筛选出该省份的数据
#         province_data = df[df[province_col] == province]
#         province_haplogroup_data = province_data[haplogroup_columns].to_numpy()
        
#         # 创建一个省份对应的哈希表，用于存储该省份中每个单倍群的计数
#         for haplogroup in haplogroups:
#             # 使用 NumPy 的向量化操作来检查单倍群是否出现在该省份的某行中
#             mask = np.any(province_haplogroup_data == haplogroup, axis=1)
            
#             # 统计出现次数（每一行最多计数一次）
#             count = np.sum(mask)
            
#             # 保存计数结果
#             haplogroup_counts[province][haplogroup] = count
    
#     return haplogroup_counts

# # 交互式输入省份字段
# province_field = input("请输入要计算的字段名称，例如 Province_cn: ")

# # 计算结果
# haplogroup_counts = calculate_haplogroup_counts(data, province_field)

# # 将结果转换为 DataFrame 以表格形式展示
# haplogroup_df = pd.DataFrame(haplogroup_counts).T.fillna(0)

# # 重置索引，使用户输入的字段成为第一列
# haplogroup_df.index.name = province_field
# haplogroup_df.reset_index(inplace=True)

# # 显示表格结果
# print(haplogroup_df)

# # 如果需要将结果保存为文件，可以使用以下命令：
# # haplogroup_df.to_csv("haplogroup_counts.csv", index=False)
