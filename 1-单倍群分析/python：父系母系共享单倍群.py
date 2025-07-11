import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

# 设置plt.rcParams以确保输出为矢量图格式并且字体为 Arial
plt.rcParams['pdf.fonttype'] = 42
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['font.family'] = 'Arial'  # 统一设置字体为Arial

# 第一步：加载第一个文件（用实际文件路径替换）
first_file_path = 'C:/Users/victo/Desktop/新建 Text Document.txt'
df = pd.read_csv(first_file_path, sep='\t', index_col=0)
df = df.apply(pd.to_numeric, errors='coerce')  # 将所有值转换为数值类型，无法转换的变为NaN

# 第二步：加载第二个文件（用实际文件路径替换）
second_file_path = 'C:/Users/victo/Desktop/新建 Text Document (2).txt'
second_df = pd.read_csv(second_file_path, sep='\t', index_col=0)
second_df = second_df.apply(pd.to_numeric, errors='coerce')

# 第三步：确保两个DataFrame（数据表）对齐，以避免索引和列错位
aligned_second_df = second_df.reindex(index=df.index, columns=df.columns)

# 第四步：创建上三角掩码和下三角掩码（包括对角线）
upper_tri_mask = np.triu(np.ones_like(df, dtype=bool))  # 上三角掩码
lower_tri_mask = np.tril(np.ones_like(aligned_second_df, dtype=bool))  # 下三角掩码

# 第五步：将两个矩阵合并
# - 使用第一个DataFrame的上三角（包括对角线）
# - 使用对齐后的第二个DataFrame的下三角（包括对角线）
combined_matrix_aligned = df.where(upper_tri_mask, aligned_second_df)

# 第六步：创建具有100种渐变颜色的自定义颜色映射
colors = ["#20364F", "#31646C", "#4E9280", "#96B89B", "#EEEFFF", "#ECD9CF", 
          "#D49C87", "#B86265", "#8B345E", "#50184E"]  # 基础颜色列表
n_colors = 100  # 颜色渐变的总数量

# 使用LinearSegmentedColormap生成颜色渐变映射
cmap = LinearSegmentedColormap.from_list("custom_cmap", colors, N=n_colors)

# 第七步：绘制组合后的热图，并调整纵横比和标签间距
plt.figure(figsize=(14, 12))  # 设置图形尺寸，使间距更大
ax = sns.heatmap(combined_matrix_aligned, cmap=cmap, annot=False, cbar=True, square=False)  # 绘制热图

# 将X轴移动到右侧
ax.xaxis.set_ticks_position('top')
ax.xaxis.set_label_position('top')

# 将Y轴保留在左侧
ax.yaxis.set_ticks_position('left')
ax.yaxis.set_label_position('left')

# 调整标签的间距和方向
plt.xticks(rotation=90, ha='left', fontsize=10)  # X轴标签
plt.yticks(rotation=0, fontsize=10)  # Y轴标签

# 设置图形标题和标签
plt.title('Sharing haplogroups of mtDNA and ChrY', fontsize=14)  # 设置图标题，字体大小14
plt.xlabel('', fontsize=12)  # X轴标签
plt.ylabel('', fontsize=12)  # Y轴标签
plt.tight_layout(pad=2.0)  # 调整布局，增加一些间距，防止标签重叠

# 显示图形
plt.show()
