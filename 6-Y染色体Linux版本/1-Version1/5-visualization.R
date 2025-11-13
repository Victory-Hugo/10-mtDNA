#!/usr/bin/env Rscript
################################################################################
# 5-visualization.R
# 数据可视化和图表生成
# 作者: BigLin
# 用法: Rscript 5-visualization.R <distance_matrices_rds> <freq_matrix_csv> <output_dir> <nc>
################################################################################

# 加载必要的库
suppressPackageStartupMessages({
    library(ggplot2, quietly=TRUE)
    library(ggrepel, quietly=TRUE)
    library(ggalluvial, quietly=TRUE)
    library(pheatmap, quietly=TRUE)
    library(ape, quietly=TRUE)
    library(scales, quietly=TRUE)
    library(dplyr, quietly=TRUE)
    library(tidyr, quietly=TRUE)
    library(grid, quietly=TRUE)
})

# 获取命令行参数
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 4) {
    stop("用法: Rscript 5-visualization.R <distance_matrices_rds> <freq_matrix_csv> <output_dir> <nc>")
}

distance_matrices_rds <- args[1]
freq_matrix_csv <- args[2]
output_dir <- args[3]
NC <- as.integer(args[4])

# 加载工具函数
source("1-utils.R")

message2("开始可视化...")
message2("距离矩阵文件: %s", distance_matrices_rds)
message2("频率矩阵文件: %s", freq_matrix_csv)
message2("输出目录: %s", output_dir)

# 读取数据
distance_matrices <- readRDS(distance_matrices_rds)
freq_mat_df <- read.csv(freq_matrix_csv, stringsAsFactors = FALSE, fileEncoding = "UTF-8")

# 创建输出目录
figs_dir <- file.path(output_dir, "figs")
dir.create(figs_dir, showWarnings = FALSE, recursive = TRUE)

# 绘图参数
plot_base_width <- 8
plot_base_height <- 6

message2("开始绘制热图...")

# 距离热图绘制函数
plot_distance_heatmap <- function(D, out_file, title="Distance heatmap", method="ward.D2") {
    cols <- nature_hm_seq(100)
    base <- tools::file_path_sans_ext(out_file)
    
    tryCatch({
        # PNG版本
        png(paste0(base, ".png"), width=plot_base_width*100, height=plot_base_height*100, res=100)
        suppressWarnings(pheatmap(D, color=cols,
                                clustering_distance_rows=as.dist(D),
                                clustering_distance_cols=as.dist(D),
                                clustering_method=method,
                                main=title))
        dev.off()
        
        # PDF版本
        pdf(paste0(base, ".pdf"), width=plot_base_width, height=plot_base_height)
        suppressWarnings(pheatmap(D, color=cols,
                                clustering_distance_rows=as.dist(D),
                                clustering_distance_cols=as.dist(D),
                                clustering_method=method,
                                main=title))
        dev.off()
        
        message2("✓ 热图保存: %s", basename(base))
    }, error = function(e) {
        message2("✗ 热图绘制失败: %s - %s", basename(base), e$message)
    })
}

# 系统发育树绘制函数
plot_dendrogram <- function(D, out_file, method="ward.D2") {
    base <- tools::file_path_sans_ext(out_file)
    
    tryCatch({
        hc <- hclust(as.dist(D), method=method)
        
        # PDF版本
        pdf(paste0(base, ".pdf"), width=plot_base_width, height=plot_base_height)
        plot(hc, xlab="", sub="", main=sprintf("Hierarchical clustering (%s)", method))
        dev.off()
        
        # PNG版本
        png(paste0(base, ".png"), width=plot_base_width*100, height=plot_base_height*100, res=100)
        plot(hc, xlab="", sub="", main=sprintf("Hierarchical clustering (%s)", method))
        dev.off()
        
        # 保存Newick格式树文件
        phy <- ape::as.phylo(hc)
        phy$tip.label <- sanitize_label_ascii(gsub_bytes("[^A-Za-z0-9_.-]", "_", phy$tip.label))
        ape::write.tree(phy, file=paste0(base, ".nwk"))
        
        message2("✓ 系统树保存: %s", basename(base))
    }, error = function(e) {
        message2("✗ 系统树绘制失败: %s - %s", basename(base), e$message)
    })
}

# 邻接树绘制函数
plot_nj_tree <- function(D, out_file, title="Neighbor-Joining tree") {
    base <- tools::file_path_sans_ext(out_file)
    
    tryCatch({
        tr <- ape::nj(as.dist(D))
        
        # PDF版本
        pdf(paste0(base, ".pdf"), width=plot_base_width, height=plot_base_height)
        plot(tr, main=title, cex=0.8)
        dev.off()
        
        # PNG版本
        png(paste0(base, ".png"), width=plot_base_width*100, height=plot_base_height*100, res=100)
        plot(tr, main=title, cex=0.8)
        dev.off()
        
        # 保存Newick格式
        tr$tip.label <- sanitize_label_ascii(gsub_bytes("[^A-Za-z0-9_.-]", "_", tr$tip.label))
        ape::write.tree(tr, file=paste0(base, ".nwk"))
        
        message2("✓ 邻接树保存: %s", basename(base))
    }, error = function(e) {
        message2("✗ 邻接树绘制失败: %s - %s", basename(base), e$message)
    })
}

# MDS分析和绘图
run_mds_plot <- function(D, out_png, title="MDS Analysis") {
    base <- tools::file_path_sans_ext(out_png)
    
    tryCatch({
        conf <- cmdscale(as.dist(D), k=2, eig=TRUE)
        pts <- as.data.frame(conf$points)
        pts$population <- sanitize_label_ascii(rownames(pts))
        colnames(pts)[1:2] <- c("Dim1", "Dim2")
        
        # 保存MDS坐标
        write.csv(pts, paste0(base, "_MDS_points.csv"), row.names=FALSE, fileEncoding="UTF-8")
        
        p <- ggplot(pts, aes(Dim1, Dim2, label=population)) +
            geom_point(size=3, color="steelblue", alpha=0.8) +
            geom_text_repel(max.overlaps = 100, size=3) +
            theme_bw(base_size=12) +
            labs(title=title, x="Dimension 1", y="Dimension 2") +
            theme(plot.title = element_text(hjust = 0.5))
        
        save_png_pdf(p, base, width=plot_base_width, height=plot_base_height)
        
        message2("✓ MDS图保存: %s", basename(base))
    }, error = function(e) {
        message2("✗ MDS绘制失败: %s - %s", basename(base), e$message)
    })
}

# Alluvial图绘制
plot_alluvial_simple_freq <- function(freq_mat_df, out_png) {
    base <- tools::file_path_sans_ext(out_png)
    
    tryCatch({
        long <- freq_mat_df %>%
            rename(population=1) %>%
            pivot_longer(-population, names_to="hap", values_to="freq")
        
        long$hap <- droplevels(factor(long$hap))
        long$population <- factor(sanitize_label_ascii(long$population))
        
        pal <- get_palette_for_levels(levels(long$hap), min_n = max(512, nlevels(long$hap)))
        
        p <- ggplot(long, aes(y=freq, axis1=hap, axis2=population)) +
            geom_alluvium(aes(fill=hap), width=1/12, knot.pos=0.4, alpha=0.9) +
            geom_stratum(aes(fill=hap), width=1/12, color="grey30") +
            geom_text(stat="stratum", aes(label=after_stat(stratum)), hjust=0, size=3) +
            scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
            scale_x_discrete(limits=c("Haplogroup","Population"), expand=c(.1,.05)) +
            scale_fill_manual(values=pal, limits=names(pal), drop=FALSE, na.value="#BDBDBD") +
            theme_bw(base_size=12) +
            labs(title = sprintf("Haplogroup Alluvial Plot (NC=%d)", NC)) +
            theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
        
        save_png_pdf(p, base, width=plot_base_width*1.6, height=plot_base_height*1.1)
        
        message2("✓ Alluvial图保存: %s", basename(base))
    }, error = function(e) {
        message2("✗ Alluvial图绘制失败: %s - %s", basename(base), e$message)
    })
}

# 主要距离类型
main_distances <- c("TVD", "JSD", "Nei_D")

# 为主要距离类型绘制所有图表
for (dist_name in main_distances) {
    D <- distance_matrices[[dist_name]]
    if (is.null(D)) next
    
    message2("绘制 %s 距离图表...", dist_name)
    
    # 热图
    plot_distance_heatmap(D, 
                         file.path(figs_dir, sprintf("heatmap_%s_NC%02d.png", dist_name, NC)),
                         title=sprintf("%s Distance Heatmap (NC=%d)", dist_name, NC))
    
    # 系统发育树
    plot_dendrogram(D, 
                   file.path(figs_dir, sprintf("dendrogram_%s_NC%02d.pdf", dist_name, NC)))
    
    # 邻接树
    plot_nj_tree(D, 
                file.path(figs_dir, sprintf("NJ_%s_NC%02d.pdf", dist_name, NC)),
                title=sprintf("%s Neighbor-Joining Tree (NC=%d)", dist_name, NC))
    
    # MDS图
    run_mds_plot(D, 
                file.path(figs_dir, sprintf("MDS_%s_NC%02d.png", dist_name, NC)),
                title=sprintf("%s MDS Analysis (NC=%d)", dist_name, NC))
}

# Alluvial图
message2("绘制 Alluvial 图...")
plot_alluvial_simple_freq(freq_mat_df, 
                         file.path(figs_dir, sprintf("alluvial_NC%02d.png", NC)))

# 频率条形图
message2("绘制频率条形图...")
tryCatch({
    long <- freq_mat_df %>%
        rename(population=1) %>%
        pivot_longer(-population, names_to="hap", values_to="freq") %>%
        arrange(desc(freq)) %>%
        group_by(population) %>%
        slice_head(n=10) %>%  # 只显示每个群体的前10个单倍群
        ungroup()
    
    pal <- get_palette_for_levels(unique(long$hap))
    
    p <- ggplot(long, aes(x=reorder(hap, freq), y=freq, fill=hap)) +
        geom_col(alpha=0.8) +
        facet_wrap(~population, scales="free_x") +
        scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
        scale_fill_manual(values=pal, na.value="#BDBDBD") +
        theme_bw(base_size=10) +
        theme(axis.text.x = element_text(angle=45, hjust=1),
              legend.position = "none",
              plot.title = element_text(hjust = 0.5)) +
        labs(x="Haplogroup", y="Frequency", 
             title=sprintf("Top Haplogroup Frequencies by Population (NC=%d)", NC))
    
    base <- file.path(figs_dir, sprintf("barplot_freq_NC%02d", NC))
    save_png_pdf(p, base, width=plot_base_width*1.5, height=plot_base_height*1.2)
    
    message2("✓ 条形图保存: barplot_freq_NC%02d", NC)
}, error = function(e) {
    message2("✗ 条形图绘制失败: %s", e$message)
})

message2("可视化完成！")
message2("输出图片数量: %d", length(list.files(figs_dir, pattern="\\.(png|pdf)$")))
message2("输出目录: %s", figs_dir)