from Bio import SeqIO
from Bio.Seq import Seq
def replace_regions_in_fasta(input_file, output_file, regions):
    """
    替换FASTA文件中指定区域的碱基为'-'
    
    :param input_file: 输入的FASTA文件路径
    :param output_file: 输出的FASTA文件路径
    :param regions: 需要替换的区域列表，每个区域用一个元组表示 (start, end)
    """
    with open(input_file, "r") as infile, open(output_file, "w") as outfile:
        for record in SeqIO.parse(infile, "fasta"):
            # 将序列转换为可修改的列表
            sequence = list(str(record.seq))
            # 替换指定区域
            for start, end in regions:
                sequence[start - 1:end] = '-' * (end - start + 1)  # start-1调整为0索引
            # 更新序列
            record.seq = Seq(''.join(sequence))
            # 写入修改后的序列
            SeqIO.write(record, outfile, "fasta")
    print(f"处理完成！修改后的文件已保存到：{output_file}")

# 设置输入文件和输出文件路径
input_file = r"C:/Users/victo/Desktop/现代DNA对齐（2024年11月26日）.fasta"
output_file = r"C:/Users/victo/Desktop/现代DNA对齐_替换后的文件.fasta"

# 定义需要替换的区域列表 [(start1, end1), (start2, end2), ...]
regions_to_replace = [
    (303, 318),  # 303到318
    (522, 531),  # 522到531
    (576, 580),   # 576到580
    (16196,16214),
    (8287,8305),
    (5907,5914),
    (3118,3120),
]

# 调用函数处理
replace_regions_in_fasta(input_file, output_file, regions_to_replace)
