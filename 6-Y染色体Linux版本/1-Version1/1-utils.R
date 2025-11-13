#!/usr/bin/env Rscript
################################################################################
# 1-utils.R
# 通用工具函数库
# 作者: BigLin
# 描述: 包含所有通用的工具函数和颜色配置
################################################################################

# 基础设置
options(stringsAsFactors = FALSE, warn = -1)
Sys.setlocale("LC_CTYPE", "en_US.UTF-8")
options(encoding = "UTF-8")

# 消息函数
message2 <- function(fmt, ...) {
    cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(fmt, ...), "\n")
}

# 单倍群清理函数
sanitize_hap <- function(x) {
    x <- stringr::str_replace_all(x, "[^A-Za-z0-9]", "")
    x <- tolower(x)
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
}

# 单倍群截断函数
hap_substr <- function(h, nc) substr(h, 1, nc)

# ASCII标签清理函数
sanitize_label_ascii <- function(x) {
    x <- as.character(x)
    y <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "_")
    y[is.na(y)] <- "_"
    y
}

# 字节安全的gsub
gsub_bytes <- function(pattern, replacement, x) {
    gsub(pattern, replacement, x, useBytes = TRUE)
}

# 颜色配置
.okabe_ito <- c("#0072B2", "#56B4E9", "#009E73", "#E69F00", "#CC79A7", "#D55E00", "#999999", "#F0E442")
.nature_seeds <- c(
    "#4477AA", "#66CCEE", "#228833", "#CCBB44", "#EE6677", "#AA3377",
    "#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77",
    "#CC6677", "#882255", "#AA4499", "#DDDDDD",
    "#1F77B4", "#FF7F0E", "#2CA02C", "#D62728", "#9467BD", "#8C564B",
    "#E377C2", "#7F7F7F", "#BCBD22", "#17BECF"
)
.seeds <- unique(c(.okabe_ito, .nature_seeds))

# 颜色转换函数
.hex_to_lab <- function(hex) {
    rgb <- t(grDevices::col2rgb(hex)/255)
    lab <- grDevices::convertColor(rgb, from="sRGB", to="Lab", scale.in=1)
    as.matrix(lab)
}

.gen_candidates <- function() {
    levL <- c(80, 72, 64, 56, 48)
    levC <- c(60, 66, 62, 58, 54)
    out <- character()
    for (i in seq_along(levL)) {
        cols <- grDevices::hcl(h = seq(0, 359, by = 1.5), c = levC[i], l = levL[i])
        out <- c(out, cols)
    }
    unique(out)
}

# Glasbey样式颜色生成
glasbey_like <- function(n, seed_colors = .seeds, max_n = 2048) {
    n <- min(n, max_n)
    seed_colors <- unique(seed_colors)
    if (n <= length(seed_colors)) return(seed_colors[seq_len(n)])
    
    cand <- setdiff(.gen_candidates(), seed_colors)
    sel <- seed_colors
    sel_lab <- .hex_to_lab(sel)
    cand_lab <- .hex_to_lab(cand)
    
    set.seed(42)
    ord <- sample.int(nrow(cand_lab))
    cand_lab <- cand_lab[ord, , drop=FALSE]
    cand <- cand[ord]
    
    min_dist_to_sel <- function(row, sel_lab) {
        min(sqrt(rowSums((t(sel_lab) - row)^2)))
    }
    
    dist_min <- apply(cand_lab, 1, min_dist_to_sel, sel_lab = sel_lab)
    
    while (length(sel) < n && nrow(cand_lab) > 0) {
        pick <- which.max(dist_min)
        sel <- c(sel, cand[pick])
        new_row <- cand_lab[pick, , drop=FALSE]
        sel_lab <- rbind(sel_lab, new_row)
        
        dnew <- sqrt(rowSums((cand_lab - matrix(new_row, nrow=nrow(cand_lab), ncol=3, byrow=TRUE))^2))
        dist_min <- pmin(dist_min, dnew)
        
        cand_lab <- cand_lab[-pick, , drop=FALSE]
        cand <- cand[-pick]
        dist_min <- dist_min[-pick]
    }
    
    sel[seq_len(n)]
}

# 为因子水平生成调色板
get_palette_for_levels <- function(fac_levels, min_n = 0) {
    lv <- as.character(unique(fac_levels))
    n <- max(length(lv), min_n)
    pal <- glasbey_like(n)
    setNames(pal[seq_len(length(lv))], lv)
}

# Nature风格热图颜色
nature_hm_seq <- function(n = 100) {
    grDevices::colorRampPalette(c(
        "#F7FBFF", "#DEEBF7", "#C6DBEF", "#9ECAE1", "#6BAED6", 
        "#4292C6", "#2171B5", "#08519C", "#08306B"
    ))(n)
}

# 保存图片函数
save_png_pdf <- function(p, file_path_noext, width, height, dpi = 400) {
    suppressWarnings({
        ggsave(paste0(file_path_noext, ".png"), p, width=width, height=height, dpi=dpi)
        ggsave(paste0(file_path_noext, ".pdf"), p, width=width, height=height)
    })
}

message2("工具函数库加载完成")