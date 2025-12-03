#!/bin/bash

# 检查是否至少提供了一个 BAM 文件路径
if [ "$#" -eq 0 ]; then
    echo "Usage: $0 <bam_file1> [bam_file2 ...]"
    exit 1
fi

# 遍历传入的 BAM 文件路径
for bam_file in "$@"; do
    # 检查文件是否存在
    if [ ! -f "$bam_file" ]; then
        echo "File $bam_file not found, skipping."
        continue
    fi

    # 从 BAM 文件名称（去掉扩展名）获取样本名称
    sample_name=$(basename "$bam_file" .bam)
    
    echo "Adding sample information to $bam_file"

    # 添加样本信息并覆盖原始文件
    samtools addreplacerg \
        -r "@RG\tID:${sample_name}\tSM:${sample_name}" \
        -o temp.bam "$bam_file" && mv temp.bam "$bam_file"

    # 为更新后的 BAM 文件创建索引
    samtools index "$bam_file"
done

echo "Sample information added and indexed for provided BAM files."
