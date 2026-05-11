import pandas as pd
import sys

# 获取传入的文件路径
input_file = sys.argv[1]
output_file = sys.argv[2]

# 定义新的颜色列表
# new_colors = [
#     '#ed3939', '#eb9237', '#e8d130',
#     '#ffff37', '#a5f130', '#4ae033',
#     '#50a161', '#56d4bb', '#3ea9b4',
#     '#308ac5', '#3d2fd8', '#a737d5','#999999'
# ]
# 定义新的颜色列表
# new_colors = [
#     '#50184E','#845ec2', '#d65db1', '#ff6f91',
#     '#ff9671', '#ffc75f', '#f9f871',
#     '#c2e577', '#6db883', '#46857a',
#     '#2481a1','#2e4b6e','#999999'
# ]
new_colors = [
    '#D55E00','#DC7900', '#E9A904', '#E0C318',
    '#f9f871', '#7AB241', '#46857a',
    '#40AECB','#238DC8','#0072B2']
# 原始颜色 #EA1F1F，饱和度降低后的颜色：#ed3939
# 原始颜色 #E88421，饱和度降低后的颜色：#eb9237
# 原始颜色 #E5C923，饱和度降低后的颜色：#e8d130
# 原始颜色 #FFF924，饱和度降低后的颜色：#ffff37
# 原始颜色 #9DEF1B，饱和度降低后的颜色：#a5f130
# 原始颜色 #42D726，饱和度降低后的颜色：#4ae033
# 原始颜色 #449657，饱和度降低后的颜色：#50a161
# 原始颜色 #4CCCB3，饱和度降低后的颜色：#56d4bb
# 原始颜色 #369BA8，饱和度降低后的颜色：#3ea9b4
# 原始颜色 #2B7EBC，饱和度降低后的颜色：#308ac5
# 原始颜色 #3626D1，饱和度降低后的颜色：#3d2fd8
# 原始颜色 #A128CE，饱和度降低后的颜色：#a737d5
# 原始颜色 #999999（灰色系列通常饱和度已很低，这里不做调整）

# 读取CSV文件
df = pd.read_csv(input_file, delimiter=';', header=None, names=['Group', 'Color', 'Value'])

# 创建一个字典来跟踪颜色的分配以避免重复
assigned_colors = {}

# 替换颜色的函数并处理重复情况
def replace_colors_v3(row, color_index):
    color = new_colors[color_index % len(new_colors)]
    if color in assigned_colors:
        occurrence = assigned_colors[color]
        row['Value'] = f'lines-{occurrence * 2 - 1}'
        assigned_colors[color] += 1
    else:
        assigned_colors[color] = 1
    row['Color'] = color
    return row

# 替换数据框中的颜色
color_index = 0
for i in range(len(df)):
    df.iloc[i] = replace_colors_v3(df.iloc[i], color_index)
    color_index += 1

# 保存修改后的数据框到新的CSV文件
df.to_csv(output_file, sep=';', index=False, header=False)
