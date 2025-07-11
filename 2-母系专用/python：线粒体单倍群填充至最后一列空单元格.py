import pandas as pd

# 读取 CSV 文件
file_path = 'C:/Users/victo/Desktop/正序等级.txt'
df = pd.read_csv(file_path,sep="\t")

# 使用 ffill 函数按行填充缺失值，填充每行最右边的空白单元格为其最右侧的非空单元格
df_filled = df.apply(lambda row: row.ffill(axis=0), axis=1)
df_filled.to_csv("C:/Users/victo/Desktop/整理.csv",index=False)


