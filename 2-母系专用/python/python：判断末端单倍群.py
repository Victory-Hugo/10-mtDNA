def find_terminal_haplogroups(input_file_path, output_file_path):
    import pandas as pd
    
    # 读取文件并忽略第一行标题
    df = pd.read_csv(input_file_path, sep='\t', header=0)
    
    # 提取所有单倍群
    all_haplogroups = set(df.iloc[:, 0])
    upstream_haplogroups = set(df.iloc[:, 1:].stack())
    
    # 找出末端单倍群（仅在第一列出现的单倍群）
    terminal_haplogroups = all_haplogroups - upstream_haplogroups
    
    # 保存结果
    with open(output_file_path, 'w') as f:
        for haplogroup in terminal_haplogroups:
            f.write(f"{haplogroup}\n")

# 输入和输出文件路径
input_file_path = 'C:/Users/victo/Desktop/逆序等级.txt'
output_file_path = 'C:/Users/victo/Desktop/terminal_haplogroups.txt'

# 调用函数处理文件
find_terminal_haplogroups(input_file_path, output_file_path)
