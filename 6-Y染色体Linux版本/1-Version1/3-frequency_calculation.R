#!/usr/bin/env Rscript
################################################################################
# 3-frequency_calculation.R
# 频率表和频率矩阵计算
# 作者: BigLin
# 用法: Rscript 3-frequency_calculation.R <input_rds> <output_dir> <nc>
################################################################################

# 加载必要的库
suppressPackageStartupMessages({
    library(dplyr, quietly=TRUE)
    library(tidyr, quietly=TRUE)
})

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("用法: Rscript 3-frequency_calculation.R <input_rds> <output_dir> <nc>")
}

input_rds <- args[1]
output_dir <- args[2]
NC <- as.integer(args[3])

# 加载工具函数
source("1-utils.R")

message2("开始频率计算...")
message2("输入数据: %s", input_rds)
message2("输出目录: %s", output_dir)

# 读取预处理后的数据
processed_data <- readRDS(input_rds)
df <- processed_data$data

message2("数据行数: %d", nrow(df))
message2("群体数: %d", length(unique(df$population)))
message2("单倍群数: %d", length(unique(df$hap)))

# 频率表计算函数
freq_table <- function(df, hap_col="hap", by_col="population") {
    df %>% 
        count(.data[[by_col]], .data[[hap_col]], name="n") %>%
        group_by(.data[[by_col]]) %>% 
        mutate(freq = n/sum(n)) %>% 
        ungroup()
}

# 频率矩阵构建函数
make_freq_matrix <- function(freq_df, by_col="population", hap_col="hap", val_col="freq") {
    freq_df %>%
        select(all_of(c(by_col, hap_col, val_col))) %>%
        tidyr::pivot_wider(
            names_from = all_of(hap_col),
            values_from = all_of(val_col),
            values_fill = 0
        ) %>%
        as.data.frame()
}

# 群体排序函数
order_pop_by_dominant_hap <- function(freq_mat_df) {
    long <- freq_mat_df %>%
        rename(population = 1) %>%
        pivot_longer(-population, names_to = "hap", values_to = "freq")
    
    pop_top <- long %>%
        group_by(population) %>%
        slice_max(freq, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        rename(top_hap = hap, top_freq = freq)
    
    hap_rank <- pop_top %>%
        group_by(top_hap) %>%
        summarise(max_freq_any_pop = max(top_freq, na.rm = TRUE), .groups = "drop") %>%
        arrange(desc(max_freq_any_pop)) %>%
        mutate(hap_order = row_number())
    
    pop_order_df <- pop_top %>%
        left_join(hap_rank, by = "top_hap") %>%
        arrange(hap_order, desc(top_freq), population)
    
    pop_order_df$population
}

# 计算频率表
message2("计算频率表...")
freq_df <- freq_table(df, hap_col="hap", by_col="population")

# 计算频率矩阵
message2("构建频率矩阵...")
freq_mat_df <- make_freq_matrix(freq_df, by_col="population", hap_col="hap", val_col="freq")

# 计算计数表和矩阵
message2("计算计数矩阵...")
count_df <- df %>% count(population, hap, name="n")
count_mat_df <- make_freq_matrix(count_df, by_col="population", hap_col="hap", val_col="n")

# 计算多样性指标
message2("计算多样性指标...")

# Shannon多样性
shannon <- function(p) { 
    p <- p[p > 0]
    -sum(p * log(p)) 
}

# Simpson多样性
simpson <- function(p) { 
    1 - sum(p^2) 
}

# 有效单倍群数
heff_shannon <- function(p) { 
    exp(shannon(p)) 
}

heff_simpson <- function(p) { 
    1/sum(p^2) 
}

# 计算每个群体的多样性
mat <- as.matrix(freq_mat_df[,-1])
rownames(mat) <- freq_mat_df$population

diversity_df <- data.frame(
    population = freq_mat_df$population,
    shannon = apply(mat, 1, shannon),
    simpson = apply(mat, 1, simpson),
    heff_shannon = apply(mat, 1, heff_shannon),
    heff_simpson = apply(mat, 1, heff_simpson),
    n_haplogroups = rowSums(mat > 0),
    total_samples = rowSums(as.matrix(count_mat_df[,-1])),
    stringsAsFactors = FALSE
)

# 群体排序
message2("计算群体排序...")
pop_order <- order_pop_by_dominant_hap(freq_mat_df)

# 创建输出目录
tables_dir <- file.path(output_dir, "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# 保存结果
message2("保存频率表和矩阵...")

write.csv(freq_df, file.path(tables_dir, sprintf("frequency_table_NC%02d.csv", NC)),
          row.names=FALSE, fileEncoding="UTF-8")

write.csv(freq_mat_df, file.path(tables_dir, sprintf("freq_matrix_NC%02d.csv", NC)),
          row.names=FALSE, fileEncoding="UTF-8")

write.csv(count_mat_df, file.path(tables_dir, sprintf("count_matrix_NC%02d.csv", NC)),
          row.names=FALSE, fileEncoding="UTF-8")

write.csv(diversity_df, file.path(tables_dir, sprintf("diversity_indices_NC%02d.csv", NC)),
          row.names=FALSE, fileEncoding="UTF-8")

# 保存群体排序
write.csv(data.frame(population = pop_order, order = seq_along(pop_order)),
          file.path(tables_dir, sprintf("population_order_NC%02d.csv", NC)),
          row.names=FALSE, fileEncoding="UTF-8")

# 创建汇总信息
summary_info <- list(
    nc = NC,
    n_populations = nrow(freq_mat_df),
    n_haplogroups = ncol(freq_mat_df) - 1,
    total_samples = sum(diversity_df$total_samples),
    processing_time = Sys.time(),
    population_order = pop_order,
    diversity_summary = summary(diversity_df[,-1]),
    top_haplogroups = names(sort(colSums(as.matrix(count_mat_df[,-1])), decreasing = TRUE))[1:min(10, ncol(count_mat_df)-1)]
)

saveRDS(summary_info, file.path(tables_dir, sprintf("frequency_summary_NC%02d.rds", NC)))

message2("频率计算完成！")
message2("输出表格数: %d", 5)
message2("群体数: %d", nrow(freq_mat_df))
message2("单倍群数: %d", ncol(freq_mat_df) - 1)
message2("总样本数: %d", sum(diversity_df$total_samples))