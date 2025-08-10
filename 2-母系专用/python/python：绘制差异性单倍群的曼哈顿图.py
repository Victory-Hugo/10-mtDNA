import pandas as pd
from collections import defaultdict
import matplotlib.pyplot as plt
from adjustText import adjust_text
import os

# 输入文件
data = pd.read_csv('C:/Users/victo/Desktop/新建 Text Document.txt', sep='\t')

# 输出文件路径
outputfile_path = "C:/Users/victo/Desktop"

# 定义函数，找到最下游的单倍群
def find_most_downstream_level(data):
    # 找到最下游的单倍群列
    level_columns = [col for col in data.columns if col.startswith('Level_')]
    level_numbers = [int(col.split('_')[1]) for col in level_columns]
    most_downstream_level = f"Level_{max(level_numbers)}"
    return most_downstream_level

# 定义函数，计算每个类别的总数
def calculate_category_totals(data, column_name):
    # 去除缺失值再计算每个类别的总数
    category_totals = data[column_name].dropna().value_counts().to_dict()
    return category_totals

# 定义函数，计算单倍群的频率
def calculate_haplogroup_frequencies(data, column_name, category_totals):
    # 动态确定最下游的单倍群
    most_downstream_level = find_most_downstream_level(data)
    
    # 清除列名中的空白并检查指定的列是否存在
    data.columns = data.columns.str.strip()
    column_name = column_name.strip()
    if (column_name not in data.columns):
        return f"列名 '{column_name}' 不存在于数据中，请检查输入。"

    # 按最下游的单倍群列和指定列进行分组，然后计算出现次数
    distribution = data.groupby([most_downstream_level, column_name]).size().unstack(fill_value=0).reset_index()

    # 获取指定列中的实际类别
    valid_categories = data[column_name].dropna().unique()  # 去除缺失值

    # 初始化一个字典来存储下游单倍群和单倍群本身
    haplogroup_totals = defaultdict(lambda: defaultdict(int))

    # 计算单倍群和下游单倍群的数量
    for haplogroup in distribution[most_downstream_level]:
        for idx, group in distribution.iterrows():
            if group[most_downstream_level].startswith(haplogroup):
                for category in valid_categories:
                    if category in group:  # 检查类别是否在行中
                        haplogroup_totals[haplogroup][category] += group.get(category, 0)

    # 转为Pandas DataFrame并保存为文件
    total_distribution = pd.DataFrame(haplogroup_totals).T.reset_index()
    total_distribution.columns = [most_downstream_level] + valid_categories.tolist()

    # 输出到桌面的txt文件
    total_distribution.to_csv(f"{outputfile_path}/Haplogroup_counts.txt", sep="\t", encoding="utf-8", index=False)

    # 计算频率
    frequencies = total_distribution.copy()
    for category in valid_categories:
        total_count = category_totals.get(category, 0)  # 使用预先计算的类别总数，并对不存在的类别返回0
        if total_count > 0:
            frequencies[category] = frequencies[category] / total_count
    return frequencies

def calculate_frequency_difference(frequency_result):
    # 询问用户输入两个列名以计算差值
    column1 = input("请输入第一个要计算差值的列名：").strip()
    column2 = input("请输入第二个要计算差值的列名：").strip()

    # 检查列是否存在于频率结果中
    if column1 not in frequency_result.columns or column2 not in frequency_result.columns:
        return f"输入的列名 '{column1}' 或 '{column2}' 不存在于频率结果中，请检查输入。"

    # 计算两个列之间的差值
    frequency_result['Difference'] = frequency_result[column1] - frequency_result[column2]
    return frequency_result, column1, column2

# 定义曼哈顿图的绘制
def plot_manhattan(frequency_result, column1, column2):
    # 确保生成的字体为矢量可编辑
    plt.rcParams['pdf.fonttype'] = 42
    plt.rcParams['ps.fonttype'] = 42

    # 提取绘图所需数据
    haplogroups = frequency_result.columns[0]
    differences = frequency_result['Difference']

    # 设定颜色和形状
    colors = ['#bb2649','#f35d74','#E88421','#E5C923','#ded82d','#8FBC8F', '#37945F',
              '#c6ffe6','#4CCCB3', '#369BA8','#2B7EBC','#d4eaf7','#003f8f','#ffc7ff','#A128CE','#6c35de','#999999','#1F1F1F' ]
    
    # 设置阈值
    thresholds = [-0.01, 0.01]
    print(f"阈值已经被设定为{thresholds}")
    plt.figure(figsize=(15, 5))

    # 为首字母分配颜色
    letter_to_color = {}
    current_color_idx = 0

    # 计算x轴位置
    x_positions = []
    for haplogroup in frequency_result[haplogroups]:
        first_letter = haplogroup[0]
        if (first_letter not in letter_to_color):
            letter_to_color[first_letter] = colors[current_color_idx % len(colors)]
            current_color_idx += 1
        x_positions.append(first_letter)

    unique_letters = list(letter_to_color.keys())
    x_ticks = [x_positions.index(letter) for letter in unique_letters]

    # 绘制单倍群的差值散点图
    texts = []
    for i, haplogroup in enumerate(frequency_result[haplogroups]):
        color = letter_to_color[haplogroup[0]]
        
        # 只绘制差异不为零的散点
        if differences[i] != 0:
            plt.scatter(i, differences[i], color=color, s=50, marker='o')
        
            # 标记超出阈值的点并附上全名标签
            if differences[i] < thresholds[0] or differences[i] > thresholds[1]:
                text = plt.text(i, differences[i], haplogroup, fontsize=8, ha='right')
                texts.append(text)
                # 画线连接点和文本标签
                plt.plot([i, i], [differences[i], differences[i]], 'k-', lw=0.5, color='gray')

    # 添加阈值线
    plt.axhline(y=thresholds[0], color='#1F1F1F', linestyle='--')
    plt.axhline(y=thresholds[1], color='#1F1F1F', linestyle='--')

    # 调整文本避免重叠
    adjust_text(texts, 
                arrowprops=dict(arrowstyle='->', color='#1F1F1F', lw=0.5, shrinkA=5),
                expand_text=(1.5, 1.5),
                expand_objects=(1.5, 1.5),
                force_text=(0.5, 0.5),
                force_objects=(0.5, 0.5),
                only_move={'points':'y', 'text':'xy'})

    # 设置x轴只显示单倍群的首字母
    plt.xticks(ticks=x_ticks, labels=unique_letters, rotation=0)
    plt.xlabel('Haplogroup')
    plt.ylabel('Frequency Difference')
    plt.title(f"Frequency Difference: {column1} - {column2}")
    plt.tight_layout()
    plt.show()

# 询问用户输入要进行分类计算的列名
column_name = input("请输入要分类计算统计的列名（例如 'Province'）：").strip()

# 计算每个类别的总数
category_totals = calculate_category_totals(data, column_name)

# 将类别总数写入一个新的txt文件
with open(f"{outputfile_path}/Category_totals.txt", 'w', encoding='utf-8') as f:
    for category, total in category_totals.items():
        f.write(f"{category}\t{total}\n")

# 计算单倍群的频率，包括所有下游单倍群
frequency_result = calculate_haplogroup_frequencies(data, column_name, category_totals)

# 显示结果并询问要计算差异的列名
if isinstance(frequency_result, pd.DataFrame):
    print(frequency_result)
    frequency_result.to_csv(f"{outputfile_path}/Haplogroup_frequency.txt", sep="\t", encoding="utf-8")
    
    # 根据用户输入计算频率差异
    frequency_difference_result, column1, column2 = calculate_frequency_difference(frequency_result)
    
    # 显示并保存频率差异结果
    if isinstance(frequency_difference_result, pd.DataFrame):
        print(frequency_difference_result)
        frequency_difference_result.to_csv(f"{outputfile_path}/Frequency_difference.txt", sep="\t", encoding="utf-8")
        
        # 绘制频率差异曼哈顿图
        plot_manhattan(frequency_difference_result, column1, column2)
    else:
        print(frequency_difference_result)
else:
    print(frequency_result)

os.remove(f"{outputfile_path}/Haplogroup_frequency.txt")
print("中间文件已经被删除了，如果你不需要删除，请你删除代码最后2行。")
