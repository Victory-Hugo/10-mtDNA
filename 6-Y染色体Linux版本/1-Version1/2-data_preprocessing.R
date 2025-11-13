#!/usr/bin/env Rscript
################################################################################
# 2-data_preprocessing.R
# 数据读取、清理和预处理
# 作者: BigLin
# 用法: Rscript 2-data_preprocessing.R <input_file> <output_rds> <nc> <use_equalize>
################################################################################

# 加载必要的库
suppressPackageStartupMessages({
    library(readxl, quietly=TRUE)
    library(dplyr, quietly=TRUE)
    library(stringr, quietly=TRUE)
})

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
    stop("用法: Rscript 2-data_preprocessing.R <input_file> <output_rds> <nc> <use_equalize> <excel_sheet>")
}

input_file <- args[1]
output_rds <- args[2]
NC <- as.integer(args[3])
use_equalize <- as.logical(args[4])
excel_sheet <- args[5]

# 加载工具函数
source("1-utils.R")

message2("开始数据预处理...")
message2("输入文件: %s", input_file)
message2("Excel工作表: %s", excel_sheet)
message2("NC值: %d", NC)
message2("等量抽样: %s", use_equalize)

# 安全读取输入文件
safe_read_input <- function(path, sheet_name) {
    ext <- tolower(tools::file_ext(path))
    
    if (ext %in% c("xlsx", "xls")) {
        df <- readxl::read_excel(path, sheet = sheet_name)
    } else if (ext == "csv") {
        df <- read.csv(path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
    } else {
        df <- read.delim(path, stringsAsFactors = FALSE, fileEncoding = "UTF-8")
    }
    
    # 列名匹配
    cn <- tolower(colnames(df))
    id_idx <- which(cn %in% c("sample_id", "sample", "id"))[1]
    hap_idx <- which(cn %in% c("haplogroup_full", "haplogroup", "hap"))[1]
    pop_idx <- which(cn %in% c("population", "pop", "group"))[1]
    
    if (any(is.na(c(id_idx, hap_idx, pop_idx)))) {
        stop("输入表必须有 sample_id / haplogroup_full / population 三列（可以大小写不同）")
    }
    
    df <- df[, c(id_idx, hap_idx, pop_idx)]
    colnames(df) <- c("sample_id", "haplogroup_full", "population")
    
    return(df)
}

# 读取原始数据
message2("读取数据文件...")
df0 <- safe_read_input(input_file, excel_sheet)

message2("原始数据: %d 行, %d 列", nrow(df0), ncol(df0))
message2("样本数: %d", nrow(df0))
message2("群体数: %d", length(unique(df0$population)))

# 数据清理
df0$hap_clean <- sanitize_hap(df0$haplogroup_full)
df0$population <- sanitize_label_ascii(df0$population)

# 移除缺失值
df0 <- df0[!is.na(df0$hap_clean) & df0$hap_clean != "" & 
           !is.na(df0$population) & df0$population != "", ]

message2("清理后数据: %d 行", nrow(df0))

# 等量抽样（如果需要）
if (use_equalize) {
    message2("进行群体等量抽样...")
    
    pop_sizes <- table(df0$population)
    target_n <- min(pop_sizes)
    
    message2("目标样本数: %d", target_n)
    
    set.seed(1)
    df_in <- do.call(rbind, lapply(split(df0, df0$population), function(tab) {
        n <- nrow(tab)
        if (n > target_n) {
            tab[sample.int(n, target_n), , drop=FALSE]
        } else {
            tab
        }
    }))
    
    message2("等量抽样后: %d 行", nrow(df_in))
} else {
    df_in <- df0
}

# 应用NC截断
message2("应用 NC=%d 截断...", NC)
df_in$hap <- hap_substr(df_in$hap_clean, NC)

# 统计信息
hap_counts <- table(df_in$hap)
pop_counts <- table(df_in$population)

message2("截断后单倍群数: %d", length(hap_counts))
message2("最常见单倍群: %s (n=%d)", names(hap_counts)[1], hap_counts[1])
message2("群体样本数分布: %s", paste(names(pop_counts), pop_counts, sep="=", collapse=", "))

# 保存处理后的数据
processed_data <- list(
    data = df_in,
    nc = NC,
    use_equalize = use_equalize,
    input_file = input_file,
    processing_time = Sys.time(),
    original_rows = nrow(df0),
    processed_rows = nrow(df_in),
    n_populations = length(unique(df_in$population)),
    n_haplogroups = length(unique(df_in$hap)),
    population_sizes = as.list(pop_counts),
    haplogroup_counts = as.list(hap_counts)
)

saveRDS(processed_data, file = output_rds)

message2("数据预处理完成！")
message2("输出文件: %s", output_rds)
message2("处理时间: %s", Sys.time())