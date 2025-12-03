#!/bin/bash

# 使用方法提示
usage() {
    echo "Usage: $0 --inputfile <input_file_path> --output_dir <output_directory> --errlorlog <error_log_file>"
    exit 1
}

# 检查参数数量是否足够（至少需要6个参数）
if [ "$#" -lt 6 ]; then
    usage
fi

# 解析命令行参数
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --inputfile)
            input_file="$2"
            shift 2
            ;;
        --output_dir)
            output_dir="$2"
            shift 2
            ;;
        --errlorlog)
            error_log="$2"
            shift 2
            ;;
        *)
            echo "未知参数: $1"
            usage
            ;;
    esac
done

# 最终命名为MT（无论是MT还是chrM都提取出来）
mt_chr="MT"

# 创建输出文件夹（如果不存在）
mkdir -p "$output_dir"

# 清空或创建错误日志文件
> "$error_log"

# 检查 samtools 是否已安装
if ! command -v samtools &> /dev/null; then
    echo "你个沙雕，samtools软件没安装啊！"
    exit 1
fi

# 检查输入的 BAM 文件是否存在且可读
if [ ! -r "$input_file" ]; then
    echo "你给的啥BAM文件?读不了啊: $input_file"
    echo "$input_file" >> "$error_log"
    exit 1
fi

# 为输入的 BAM 文件生成索引（如果索引文件不存在）
if [ ! -f "${input_file}.bai" ]; then
    echo "Indexing BAM file: $input_file"
    samtools index "$input_file"
    if [ $? -ne 0 ]; then
        echo "建立索引失败,这种问题都让你碰上了: $input_file"
        echo "$input_file" >> "$error_log"
        exit 1
    fi
else
    echo "BAM index already exists for: $input_file"
fi

# 提取线粒体DNA（同时考虑MT和chrM两种情况），并最终命名为MT
output_file="$output_dir/$(basename ${input_file%.bam}_${mt_chr}.bam)"
echo "Extracting mitochondrial DNA (MT and chrM) to: $output_file"
samtools view -b "$input_file" "$mt_chr" "chrM" > "$output_file"
if [ $? -ne 0 ]; then
    echo "提取mtDNA失败了,我也不知道为什么 $input_file"
    echo "$input_file" >> "$error_log"
    exit 1
fi

# 为生成的线粒体 BAM 文件建立索引
echo "Indexing mitochondrial BAM file: $output_file"
samtools index "$output_file"
if [ $? -ne 0 ]; then
    echo "建立索引失败,这种问题都让你碰上了: $output_file"
    echo "$input_file" >> "$error_log"
    exit 1
fi

# 打印完成信息
echo "任务完成!输出文件在: $output_file. 如果有错误,都放在 $error_log."
