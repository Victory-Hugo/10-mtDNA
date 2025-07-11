import pandas as pd

# 读取参考序列
reference_file_path = 'C:/Users/victo/Desktop/GRCH38线粒体DNA.fasta'
with open(reference_file_path, 'r') as ref_file:
    reference_seq = ''.join(line.strip() for line in ref_file if not line.startswith('>'))

# 读取变异数据
data_file_path = 'C:/Users/victo/Desktop/MBE-18-0597.R2-Supplementary-Tables.csv'
data = pd.read_csv(data_file_path, skiprows=1)

# 处理所有样本
output_file_path = 'C:/Users/victo/Desktop/restored_samples.fasta'

with open(output_file_path, 'w') as fasta_file:
    for index, row in data.iterrows():
        sample_id = row['Sample ID']
        variants = row['Variants']
        
        # 创建样本序列
        sample_seq = list(reference_seq)
        
        if pd.notna(variants):
            for variant in variants.split():
                pos = int(''.join(filter(str.isdigit, variant[:-1]))) - 1  # 转换为0-based索引
                base = variant[-1]
                sample_seq[pos] = base
        
        sample_seq_str = ''.join(sample_seq)
        
        # 写入fasta文件
        fasta_file.write(f'>{sample_id}\n')
        fasta_file.write(f'{sample_seq_str}\n')

print(f"还原的fasta文件已保存为 {output_file_path}")
