#!/bin/bash
# 基础路径设置
FASTA_PATH="/mnt/e/Scientifc_software/fastHaN-main/data/merged_biallelic_AGCT_SNP_Full.fasta"
OUTPUT_PATH="/mnt/e/Scientifc_software/fastHaN-main/output"
SOFTWARE_PATH="/mnt/e/Scientifc_software/fastHaN-main/fastHaN-main/"
META_FILE="/mnt/e/Scientifc_software/fastHaN-main/data/M7b1a单倍群_META.tsv"

# 获取META文件的列数（不包括第一列）
META_COLS=$(awk '{print NF-1; exit}' "$META_FILE")

# 动态生成META_DIRS数组
META_DIRS=()
for ((i=1; i<=META_COLS; i++)); do
    META_DIRS+=("Column_$i")
done

# 创建输出文件夹
mkdir -p "$OUTPUT_PATH"
for dir in "ORIGIN_TCS" "MJN" "TCS" "MSN" "${META_DIRS[@]}"; do
    mkdir -p "$OUTPUT_PATH/$dir"
done

# 拆分META文件
split_meta_files() {
    input_file=$1
    output_dir=$2
    # 读取每一行（除去第一行标题）
    while IFS= read -r line || [[ -n "$line" ]]; do
        id=$(echo "$line" | awk '{print $1}')
        # 获取列数（包括第一列）
        num_columns=$(awk '{print NF}' <<< "$line")
        # 遍历每一列
        for ((i=2; i<=num_columns; i++)); do
            value=$(echo "$line" | awk -v col="$i" '{print $col}')
            echo -e "$id\t$value" >> "${output_dir}/Column_$((i-1)).txt"
        done
    done < <(tail -n +2 "$input_file")
}

split_meta_files "$META_FILE" "$OUTPUT_PATH"

# 输出的 FASTA 文件路径
OUTPUT_FASTA="${OUTPUT_PATH}/filtered_sequences.fasta"

# 提取序列并保存
python3 "${SOFTWARE_PATH}extract_fasta_sequences.py" "$META_FILE" "$FASTA_PATH" "$OUTPUT_FASTA"

# 文件转换并压缩
GZ_FILE="$OUTPUT_PATH/fasta2phylip.phy.gz"
python3 "${SOFTWARE_PATH}Pipline/Fasta2Phylip.py" "$OUTPUT_FASTA" | gzip > "$GZ_FILE"

# 定义一个函数来执行网络计算
run_network_algorithm() {
    algorithm=$1
    input_file=$2
    output_dir=$3
    extra_args=$4
    "${SOFTWARE_PATH}fastHaN_linux" $algorithm -i "$input_file" $extra_args -o "$output_dir"
}

# 执行各种网络计算算法
run_network_algorithm "original_tcs" "$GZ_FILE" "$OUTPUT_PATH/ORIGIN_TCS" "-t 16 -a 1 -m 1"
run_network_algorithm "mjn" "$GZ_FILE" "$OUTPUT_PATH/MJN" "-t 16 -e 0"
run_network_algorithm "modified_tcs" "$GZ_FILE" "$OUTPUT_PATH/TCS" "-t 16"
run_network_algorithm "msn" "$GZ_FILE" "$OUTPUT_PATH/MSN" "-e 0"

# 生成网络配置文件
generate_network_config() {
    json_file=$1
    meta_file=$2
    output_dir=$3
    python3 "${SOFTWARE_PATH}Script/GenNetworkConfig.py" "$json_file" "$meta_file" "$output_dir"
}

for algorithm in "ORIGIN_TCS" "MJN" "TCS" "MSN"; do
    for dir in "${META_DIRS[@]}"; do
        generate_network_config "$OUTPUT_PATH/${algorithm}.json" "$OUTPUT_PATH/${dir}.txt" "$OUTPUT_PATH/$dir"
    done
done

# 替换颜色配置文件
replace_color_config() {
    input_file=$1
    output_file=$2
    python3 "${SOFTWARE_PATH}color_replace.py" "$input_file" "$output_file"
    echo "正在为你的群组或单倍群创建彩虹色"
}

for dir in "${META_DIRS[@]}"; do
    replace_color_config "$OUTPUT_PATH/${dir}_groupconf.csv" "$OUTPUT_PATH/${dir}_groupconf_updated.csv"
done
