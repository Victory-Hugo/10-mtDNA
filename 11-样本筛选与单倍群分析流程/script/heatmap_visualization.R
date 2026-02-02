#!/usr/bin/env Rscript
# ============================================================================
# 频率矩阵热图可视化
# 使用 pheatmap 绘制聚类热图，支持分组标注
# 参数化版本：接收命令行参数以支持管道集成
# 
# 使用说明：包含多个绘图版本，使用者可自行选择需要的版本进行取消注释
# ============================================================================

# 加载必要的包
library(ape)
library(igraph)
library(ggplot2)
library(pheatmap)
library(reshape2)
library(ggsci)
library(gridExtra)

# ============================================================================
# 1. 命令行参数解析
# ============================================================================

args <- commandArgs(trailingOnly = TRUE)

# 定义参数解析函数
parse_args <- function(args) {
  params <- list()
  i <- 1
  while (i <= length(args)) {
    if (grepl("^--", args[i])) {
      key <- sub("^--", "", args[i])
      # 检查是否为 --key=value 格式
      if (grepl("=", key)) {
        kv <- strsplit(key, "=", fixed = TRUE)[[1]]
        params[[kv[1]]] <- kv[2]
      } else if (i + 1 <= length(args)) {
        # 检查下一个参数是否是值（不以 -- 开头）
        params[[key]] <- args[i + 1]
        i <- i + 1
      }
    }
    i <- i + 1
  }
  return(params)
}

params <- parse_args(args)

# 规范化参数名（将 dash 转为 underscore）
params_normalized <- list()
for (name in names(params)) {
  normalized_name <- gsub("-", "_", name)
  params_normalized[[normalized_name]] <- params[[name]]
}
params <- params_normalized

# 验证必需参数
if (is.null(params$frequency_matrix) || is.null(params$group_mapping) || is.null(params$output_dir)) {
  cat("使用方法:\n")
  cat("Rscript heatmap_visualization.R \\\n")
  cat("  --frequency-matrix <frequency_csv_path> \\\n")
  cat("  --group-mapping <group_mapping_csv_path> \\\n")
  cat("  --output-dir <output_directory> \\\n")
  cat("  --output-prefix <output_file_prefix> (可选，默认: heatmap)\n")
  quit(status = 1)
}

# 设置默认输出前缀
output_prefix <- params$output_prefix
if (is.null(output_prefix)) {
  output_prefix <- "heatmap"
}

# 创建输出目录
if (!dir.exists(params$output_dir)) {
  dir.create(params$output_dir, showWarnings = FALSE, recursive = TRUE)
}

# 构建文件路径
frequency_file <- params$frequency_matrix
group_file <- params$group_mapping
output_dir <- params$output_dir

cat("频率矩阵文件:", frequency_file, "\n")
cat("分组映射文件:", group_file, "\n")
cat("输出目录:", output_dir, "\n")
cat("输出文件前缀:", output_prefix, "\n")

# ============================================================================
# 2. 数据读取
# ============================================================================

cat("正在读取数据...\n")

# 读取频率矩阵数据
mydata <- read.table(frequency_file, header = TRUE, sep = ",", row.names = 1)
cat("频率矩阵维度:", nrow(mydata), "x", ncol(mydata), "\n")

# 读取分组标注
group <- read.table(group_file, header = TRUE, sep = ",", row.names = 1)
cat("分组标注维度:", nrow(group), "x", ncol(group), "\n")

# ============================================================================
# 3. 热图绘制 - 多个版本
# ============================================================================

cat("开始绘制热图...\n")

# ============================================================================
# 版本1：简洁风格 - 无聚类
# ============================================================================
pdf(file.path(output_dir, paste0(output_prefix, "_v1_simple.pdf")), width = 11, height = 8.5)

pheatmap(mydata,
         cluster_cols = FALSE,           # 不进行列聚类
         cluster_rows = FALSE,           # 不进行行聚类
         angle_col = c("90"),            # X轴文字角度
         fontsize = 1,                   # 字体大小
         fontsize_row = 1,               # 行标签字体大小
         fontsize_col = 1,               # 列标签字体大小
         annotation_row = group,         # 行标注
         cellwidth = 1,                  # 单元格宽度
         cellheight = 1,                 # 单元格高度
         cutree_cols = 8,                # 列聚类数量
         cutree_rows = 4,                # 行聚类数量
         main = "z-score matrix",
         color = colorRampPalette(c("#20364F", "#31646C", "#4E9280", "#96B89B", "#DCDFD2", 
                                     "#ECD9CF", "#D49C87", "#B86265", "#8B345E", "#50184E"))(100),
         lwd = 0.10)

dev.off()

cat("✅ 版本1 (简洁) 已保存\n")

# # ============================================================================
# # 版本2：多彩聚类 - 带符号标记
# # ============================================================================
# pdf(file.path(output_dir, paste0(output_prefix, "_v2_clustered_multicolor.pdf")), width = 11, height = 8.5)

# pheatmap(mydata,
#          cluster_cols = TRUE,
#          cluster_rows = TRUE,
#          angle_col = c("45"),
#          fontsize = 8,
#          fontsize_row = 8,
#          fontsize_col = 6,
#          annotation_col = group,
#          annotation_row = group,
#          cellwidth = 8,
#          cellheight = 8,
#          cutree_cols = 4,
#          cutree_rows = 4,
#          main = "FstMatrix",
#          color = colorRampPalette(c("#023047", "#126883", "#279EBC", "#90C9E6", "#FC9E7F", 
#                                      "#F75B41", "#D52120", "#20364F", "#31646C", "#4E9280", 
#                                      "#96B89B"))(10000),
#          display_numbers = matrix(ifelse(abs(mydata) > 50, "++", 
#                                          ifelse(abs(mydata) >= 40, "+", " ")), nrow(mydata)))

# dev.off()

# cat("✅ 版本2 (多彩聚类) 已保存\n")

# # ============================================================================
# # 版本3：蓝红配色 - 聚类版
# # ============================================================================
# pdf(file.path(output_dir, paste0(output_prefix, "_v3_blueRed_clustered.pdf")), width = 11, height = 8.5)

# pheatmap(mydata,
#          cluster_cols = TRUE,
#          cluster_rows = TRUE,
#          angle_col = c("90"),
#          fontsize = 8,
#          fontsize_row = 8,
#          fontsize_col = 8,
#          annotation_col = group,
#          annotation_row = group,
#          cellwidth = 8,
#          cellheight = 8,
#          cutree_cols = 4,
#          cutree_rows = 4,
#          main = "FstMatrix",
#          color = colorRampPalette(c("#023047", "#DCDFD2", "#B86265"))(10000),
#          display_numbers = matrix(ifelse(abs(mydata) > 50, "++", 
#                                          ifelse(abs(mydata) >= 40, "+", " ")), nrow(mydata)))

# dev.off()

# cat("✅ 版本3 (蓝红聚类) 已保存\n")

# # ============================================================================
# # 版本4：绿蓝配色 - 聚类版
# # ============================================================================
# pdf(file.path(output_dir, paste0(output_prefix, "_v4_greenBlue_clustered.pdf")), width = 11, height = 8.5)

# pheatmap(mydata,
#          cluster_cols = TRUE,
#          cluster_rows = TRUE,
#          angle_col = c("45"),
#          fontsize = 8,
#          fontsize_row = 8,
#          fontsize_col = 6,
#          annotation_col = group,
#          annotation_row = group,
#          cellwidth = 8,
#          cellheight = 8,
#          cutree_cols = 4,
#          cutree_rows = 4,
#          main = "FstMatrix",
#          color = colorRampPalette(c("#F8F8FF", "#91D1C2FF", "#3C5488FF"))(10000),
#          display_numbers = matrix(ifelse(abs(mydata) > 50, "++", 
#                                          ifelse(abs(mydata) >= 40, "+", " ")), nrow(mydata)))

# dev.off()

# cat("✅ 版本4 (绿蓝聚类) 已保存\n")

# ============================================================================
# 执行完成
# ============================================================================

cat("\n✅ 所有热图版本已生成完成！\n")