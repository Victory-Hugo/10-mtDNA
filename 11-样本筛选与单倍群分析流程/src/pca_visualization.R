

# ========================== 参数解析 =============================
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
    idx <- which(args == flag)
    if (length(idx) == 1 && idx < length(args)) {
        return(args[idx + 1])
    }
    return(default)
}

PCA_result_file <- get_arg("--pca-input")
color_file <- get_arg("--color-file")
output_dir <- get_arg("--output-dir", ".")
output_name <- get_arg("--output-name", "PCA_visualization.pdf")
if (is.null(PCA_result_file) || is.null(color_file)) {
    cat("[参数错误] --pca-input 和 --color-file 必须指定\n")
    quit(status = 1)
}
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
OUT_pdf_file <- file.path(output_dir, output_name)

# ========================== 主流程 =============================
library(tidyverse)
library(RColorBrewer)

cat("正在读取PCA数据...\n")
data <- read.csv(PCA_result_file, header = TRUE)

cat("数据概览:\n")
cat("样本数量:", nrow(data), "\n")
cat("变量数量:", ncol(data), "\n")
print(str(data))

# 提取分类信息
Class_big_levels <- unique(data$Class_big)
class_small_levels <- unique(data$Class_small)

cat("\nClass_big类别 (", length(Class_big_levels), "个):", paste(Class_big_levels, collapse=", "), "\n")
cat("Class_small类别 (", length(class_small_levels), "个):", paste(class_small_levels, collapse=", "), "\n")

# 提取主成分数据
PC1 <- data$PC1
PC2 <- data$PC2

frame <- data.frame(
    ID = data$ID,
    PC1 = PC1,
    PC2 = PC2,
    Class_big = data$Class_big,
    Class_small = data$Class_small
)

cat("\n数据框创建完成，包含", nrow(frame), "个样本\n")

# 读取颜色映射
color_df <- read.csv(color_file, header = TRUE, stringsAsFactors = FALSE)
color_map_from_file <- setNames(color_df$color, color_df$Class_big)
auto_palette <- brewer.pal(max(8, length(Class_big_levels)), "Set2")
if (length(Class_big_levels) > length(auto_palette)) {
    auto_palette <- colorRampPalette(brewer.pal(8, "Set2"))(length(Class_big_levels))
}
auto_palette_idx <- 1
color_mapping <- c()
for (big_class in Class_big_levels) {
    if (!is.na(color_map_from_file[big_class])) {
        color_mapping[big_class] <- color_map_from_file[big_class]
    } else {
        color_mapping[big_class] <- auto_palette[auto_palette_idx]
        auto_palette_idx <- auto_palette_idx + 1
        cat(sprintf("%s类别没有分配颜色，请注意\n", big_class))
    }
}

shape_palette <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25)
create_shape_mapping <- function(data, big_levels, shape_values) {
    shape_map <- c()
    for (big_class in big_levels) {
        small_classes <- unique(data$Class_small[data$Class_big == big_class])
        small_classes <- sort(small_classes)
        shape_idx <- 1
        for (small_class in small_classes) {
            shape_map[small_class] <- shape_values[shape_idx]
            shape_idx <- shape_idx + 1
            if (shape_idx > length(shape_values)) shape_idx <- 1
        }
    }
    return(shape_map)
}
shape_mapping <- create_shape_mapping(data, Class_big_levels, shape_palette)

cat("\n=== 颜色映射 ===\n")
for (i in seq_along(color_mapping)) {
    cat(sprintf("%-12s: %s\n", names(color_mapping)[i], color_mapping[i]))
}
cat("\n=== 形状映射 ===\n")
for (i in seq_along(shape_mapping)) {
    cat(sprintf("%-20s: %2d\n", names(shape_mapping)[i], shape_mapping[i]))
}

plot_params <- list(
    pdf_width = 35,
    pdf_height = 35,
    axis_text_size = 2,
    axis_label_size = 3,
    title_size = 3,
    point_size = 3,
    point_alpha = 0.8,
    point_lwd = 3,
    legend_text_size = 2,
    legend_point_size = 2,
    legend_ncol = 5
)

cat("\n=== 绘图参数 ===\n")
cat("散点大小:", plot_params$point_size, "\n")
cat("坐标轴文字大小:", plot_params$axis_text_size, "\n")
cat("坐标轴标签大小:", plot_params$axis_label_size, "\n")

cat("\n开始绘制PCA散点图...\n")
pdf(OUT_pdf_file, width = plot_params$pdf_width, height = plot_params$pdf_height)
layout(matrix(c(1, 2)), widths = c(1, 1), heights = c(4, 1))
par(mar = c(5, 5, 4, 2))
plot(PC1, PC2, 
     type = "n",                           # 不绘制点，只创建框架
     xlab = "PC1",                         # X轴标签
     ylab = "PC2",                         # Y轴标签
     main = "PCA",           # 标题
     cex.axis = plot_params$axis_text_size,    # 坐标轴刻度文字大小
     cex.lab = plot_params$axis_label_size,    # 坐标轴标签大小
     cex.main = plot_params$title_size,        # 标题大小
     mgp = c(3, 1, 0))                     # 坐标轴位置调整

# 按分类绘制散点，支持描边粗细参数
draw_points_by_class <- function(frame, big_levels, color_map, shape_map, point_size, point_lwd = 1) {
    points_drawn <- 0
    for (big_class in big_levels) {
        small_classes <- unique(frame$Class_small[frame$Class_big == big_class])
        small_classes <- sort(small_classes)
        cat("绘制", big_class, "类别，包含", length(small_classes), "个子类别\n")
        for (small_class in small_classes) {
            current_data <- subset(frame, Class_big == big_class & Class_small == small_class)
            if (nrow(current_data) > 0) {
                points(current_data$PC1, current_data$PC2,
                       pch = shape_map[small_class],
                       bg = color_map[big_class],
                       col = color_map[big_class],
                       cex = point_size,
                       lwd = point_lwd)
                points_drawn <- points_drawn + nrow(current_data)
                cat("  -", small_class, ":", nrow(current_data), "个点\n")
            }
        }
    }
    cat("总共绘制了", points_drawn, "个散点\n")
}

draw_points_by_class(frame, Class_big_levels, color_mapping, shape_mapping, plot_params$point_size, plot_params$point_lwd)

# === 绘制图例 ===
par(mar = c(1, 1, 1, 1))
plot.new()
create_legend_data <- function(data, big_levels, color_map, shape_map) {
    legend_items <- list(
        labels = c(),
        colors = c(),
        shapes = c(),
        fonts = c()
    )
    for (big_class in big_levels) {
        legend_items$labels <- c(legend_items$labels, big_class)
        legend_items$colors <- c(legend_items$colors, NA)
        legend_items$shapes <- c(legend_items$shapes, NA)
        legend_items$fonts <- c(legend_items$fonts, 2)
        small_classes <- unique(data$Class_small[data$Class_big == big_class])
        small_classes <- sort(small_classes)
        for (small_class in small_classes) {
            legend_items$labels <- c(legend_items$labels, paste("  ", small_class))
            legend_items$colors <- c(legend_items$colors, color_map[big_class])
            legend_items$shapes <- c(legend_items$shapes, shape_map[small_class])
            legend_items$fonts <- c(legend_items$fonts, 1)
        }
    }
    return(legend_items)
}

legend_data <- create_legend_data(data, Class_big_levels, color_mapping, shape_mapping)

legend("center", 
       legend = legend_data$labels,
       pch = legend_data$shapes,
       pt.bg = legend_data$colors,
       col = legend_data$colors,
       ncol = plot_params$legend_ncol,
       cex = plot_params$legend_text_size,
       pt.cex = plot_params$legend_point_size,
       bty = "n",
       text.font = legend_data$fonts,
       xpd = TRUE)


dev.off()
cat("\nPDF图形已保存为:", OUT_pdf_file, "\n")
quit(status = 0)
