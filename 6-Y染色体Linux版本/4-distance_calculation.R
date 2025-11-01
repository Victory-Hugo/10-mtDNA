#!/usr/bin/env Rscript
################################################################################
# 4-distance_calculation.R
# 群体间距离矩阵计算
# 作者: BigLin
# 用法: Rscript 4-distance_calculation.R <freq_matrix_csv> <output_dir> <nc>
################################################################################

# 加载必要的库
suppressPackageStartupMessages({
    library(vegan, quietly=TRUE)
    library(ape, quietly=TRUE)
})

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 3) {
    stop("用法: Rscript 4-distance_calculation.R <freq_matrix_csv> <output_dir> <nc>")
}

freq_matrix_csv <- args[1]
output_dir <- args[2]
NC <- as.integer(args[3])

# 加载工具函数
source("1-utils.R")

message2("开始距离计算...")
message2("频率矩阵文件: %s", freq_matrix_csv)
message2("输出目录: %s", output_dir)

# 读取频率矩阵
freq_mat_df <- read.csv(freq_matrix_csv, stringsAsFactors = FALSE, fileEncoding = "UTF-8")

message2("群体数: %d", nrow(freq_mat_df))
message2("单倍群数: %d", ncol(freq_mat_df) - 1)

# 转换为矩阵
mat <- as.matrix(freq_mat_df[,-1])
rownames(mat) <- freq_mat_df[,1]

message2("矩阵维度: %d x %d", nrow(mat), ncol(mat))

# 距离计算函数
# Total Variation Distance
TVD <- function(p, q) { 
    0.5 * sum(abs(p - q)) 
}

# Jensen-Shannon Distance
JSD_pair <- function(p, q, eps = 1e-12) {
    p <- p + eps
    p <- p / sum(p)
    q <- q + eps
    q <- q / sum(q)
    m <- 0.5 * (p + q)
    
    KL <- function(a, b) sum(a * log(a / b))
    jsd <- 0.5 * KL(p, m) + 0.5 * KL(q, m)
    sqrt(jsd)
}

# Nei's genetic distance
.nei_I <- function(p, q) {
    num <- sum(p * q)
    den <- sqrt(sum(p^2) * sum(q^2))
    if (den == 0) return(NA_real_)
    num / den
}

.nei_D <- function(p, q) {
    I <- .nei_I(p, q)
    if (is.na(I) || I <= 0) return(NA_real_)
    -log(I)
}

# Cavalli-Sforza和Edwards弦距离
.DA_pair <- function(p, q) { 
    1 - sum(sqrt(p * q)) 
}

# 成对距离矩阵计算
pairwise_matrix <- function(mat, FUN) {
    n <- nrow(mat)
    D <- matrix(0, n, n)
    
    for (i in seq_len(n)) {
        for (j in seq_len(n)) {
            D[i, j] <- FUN(mat[i, ], mat[j, ])
        }
    }
    
    rownames(D) <- rownames(mat)
    colnames(D) <- rownames(mat)
    D
}

# 计算各种距离
message2("计算 TVD 距离矩阵...")
D_TVD <- pairwise_matrix(mat, TVD)

message2("计算 Jensen-Shannon 距离矩阵...")
D_JSD <- pairwise_matrix(mat, function(a,b) JSD_pair(a, b))

message2("计算 Nei's 遗传距离矩阵...")
D_Nei <- pairwise_matrix(mat, .nei_D)

message2("计算 Cavalli-Sforza 和 Edwards 弦距离矩阵...")
D_DA <- pairwise_matrix(mat, .DA_pair)

# 计算欧氏距离和曼哈顿距离
message2("计算标准距离矩阵...")
D_Euclidean <- as.matrix(dist(mat, method = "euclidean"))
D_Manhattan <- as.matrix(dist(mat, method = "manhattan"))

# 计算Bray-Curtis距离（使用vegan包）
D_BrayCurtis <- as.matrix(vegdist(mat, method = "bray"))

# 创建输出目录
tables_dir <- file.path(output_dir, "tables")
dir.create(tables_dir, showWarnings = FALSE, recursive = TRUE)

# 保存距离矩阵
message2("保存距离矩阵...")

write.csv(D_TVD, file.path(tables_dir, sprintf("pairwise_TVD_NC%02d.csv", NC)),
          fileEncoding="UTF-8")

write.csv(D_JSD, file.path(tables_dir, sprintf("pairwise_JSDist_NC%02d.csv", NC)),
          fileEncoding="UTF-8")

write.csv(D_Nei, file.path(tables_dir, sprintf("pairwise_NeiD_NC%02d.csv", NC)),
          fileEncoding="UTF-8")

write.csv(D_DA, file.path(tables_dir, sprintf("pairwise_DA_NC%02d.csv", NC)),
          fileEncoding="UTF-8")

write.csv(D_Euclidean, file.path(tables_dir, sprintf("pairwise_Euclidean_NC%02d.csv", NC)),
          fileEncoding="UTF-8")

write.csv(D_Manhattan, file.path(tables_dir, sprintf("pairwise_Manhattan_NC%02d.csv", NC)),
          fileEncoding="UTF-8")

write.csv(D_BrayCurtis, file.path(tables_dir, sprintf("pairwise_BrayCurtis_NC%02d.csv", NC)),
          fileEncoding="UTF-8")

# 计算距离矩阵的统计汇总
distance_summary <- data.frame(
    Distance_Method = c("TVD", "Jensen_Shannon", "Nei_D", "DA", "Euclidean", "Manhattan", "Bray_Curtis"),
    Mean = c(
        mean(D_TVD[upper.tri(D_TVD)], na.rm = TRUE),
        mean(D_JSD[upper.tri(D_JSD)], na.rm = TRUE),
        mean(D_Nei[upper.tri(D_Nei)], na.rm = TRUE),
        mean(D_DA[upper.tri(D_DA)], na.rm = TRUE),
        mean(D_Euclidean[upper.tri(D_Euclidean)], na.rm = TRUE),
        mean(D_Manhattan[upper.tri(D_Manhattan)], na.rm = TRUE),
        mean(D_BrayCurtis[upper.tri(D_BrayCurtis)], na.rm = TRUE)
    ),
    Median = c(
        median(D_TVD[upper.tri(D_TVD)], na.rm = TRUE),
        median(D_JSD[upper.tri(D_JSD)], na.rm = TRUE),
        median(D_Nei[upper.tri(D_Nei)], na.rm = TRUE),
        median(D_DA[upper.tri(D_DA)], na.rm = TRUE),
        median(D_Euclidean[upper.tri(D_Euclidean)], na.rm = TRUE),
        median(D_Manhattan[upper.tri(D_Manhattan)], na.rm = TRUE),
        median(D_BrayCurtis[upper.tri(D_BrayCurtis)], na.rm = TRUE)
    ),
    Max = c(
        max(D_TVD[upper.tri(D_TVD)], na.rm = TRUE),
        max(D_JSD[upper.tri(D_JSD)], na.rm = TRUE),
        max(D_Nei[upper.tri(D_Nei)], na.rm = TRUE),
        max(D_DA[upper.tri(D_DA)], na.rm = TRUE),
        max(D_Euclidean[upper.tri(D_Euclidean)], na.rm = TRUE),
        max(D_Manhattan[upper.tri(D_Manhattan)], na.rm = TRUE),
        max(D_BrayCurtis[upper.tri(D_BrayCurtis)], na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
)

write.csv(distance_summary, file.path(tables_dir, sprintf("distance_summary_NC%02d.csv", NC)),
          row.names=FALSE, fileEncoding="UTF-8")

# 保存距离矩阵列表用于后续分析
distance_matrices <- list(
    TVD = D_TVD,
    JSD = D_JSD,
    Nei_D = D_Nei,
    DA = D_DA,
    Euclidean = D_Euclidean,
    Manhattan = D_Manhattan,
    BrayCurtis = D_BrayCurtis
)

saveRDS(distance_matrices, file.path(tables_dir, sprintf("distance_matrices_NC%02d.rds", NC)))

message2("距离计算完成！")
message2("计算的距离类型数: %d", nrow(distance_summary))
message2("群体对数: %d", sum(upper.tri(D_TVD)))
message2("距离汇总文件已保存")