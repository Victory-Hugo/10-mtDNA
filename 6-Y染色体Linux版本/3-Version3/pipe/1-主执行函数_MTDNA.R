###############################################
###############################################
# YHCtool_patched_sig.R
# One-click population genetics analysis (final HH-based, v2, patched + significance)
# 
# What’s new in this patched build:
# 1) Multilayer alluvial fixes:
#    - Keep all axis variables as character to prevent numeric stratum labels.
#    - Legend shows only the deepest haplogroup level; population legend removed.
#    - Text labeling: population names always on the rightmost axis; other hap axes
#      show labels only when the mean within-population frequency >= label_min_freq (default 0.01).
# 2) Unified quasirandom jitter for three scatter-like plots: method="pseudorandom", width=0.22, alpha=0.60, size=0.9.
# 3) Bootstrap module:
#    - Optional Wilcoxon pairwise significance annotations on hap-frequency plots (per hap, across populations).
#    - BH multiple-testing correction; user controls for showing significance, alpha, max brackets, and a reference population.
# 4) All ggsave calls set limitsize = FALSE to avoid the 50-inch size error for large facet or sankey graphics.
#
# Usage (minimal):
#   Rscript YHCtool_patched_sig.R
# Required inputs (same directory unless absolute paths given):
#   - input_xlsx: Excel with columns [sample_id, haplogroup_full, population] (order flexible; auto-detected)
# Optional inputs:
#   - meta_csv: population metadata with 'population' and 'group' or 'region' (for coloring/grouping, geo tests)
#   - ystr_csv: Y-STR table including numeric loci columns to compute pairwise Rst across populations
#
# Export root directory: root_out (default "PPG_YHyploResults")
# Contact: set bug reports to your internal channel
###############################################
###############################################

options(stringsAsFactors = FALSE)

# --- Parse command-line arguments ---
args <- commandArgs(trailingOnly = TRUE)
wd_path <- NULL
input_xlsx <- NULL
meta_csv <- NULL
ystr_csv <- NULL
haplogroup_hierarchy_csv <- NULL
root_out <- NULL

# Parse arguments in key=value format
for (arg in args) {
  if (grepl("^wd=", arg)) {
    wd_path <- sub("^wd=", "", arg)
  } else if (grepl("^input_xlsx=", arg)) {
    input_xlsx <- sub("^input_xlsx=", "", arg)
  } else if (grepl("^meta_csv=", arg)) {
    meta_csv <- sub("^meta_csv=", "", arg)
  } else if (grepl("^ystr_csv=", arg)) {
    ystr_csv <- sub("^ystr_csv=", "", arg)
  } else if (grepl("^haplogroup_hierarchy_csv=", arg)) {
    haplogroup_hierarchy_csv <- sub("^haplogroup_hierarchy_csv=", "", arg)
  } else if (grepl("^root_out=", arg)) {
    root_out <- sub("^root_out=", "", arg)
  }
}

# --- Working directory ---
if (!is.null(wd_path)) {
  setwd(wd_path)
} else {
  setwd("/mnt/f/OneDrive/文档（科研）/脚本/Download/10-mtDNA/6-Y染色体Linux版本/3-Version3/pipe")
}

# --- Input files (with defaults) ---
if (is.null(input_xlsx)) input_xlsx <- "../input/Table_Haplogroup.xlsx"
if (is.null(meta_csv))   meta_csv   <- "../conf/Group.csv"
if (is.null(ystr_csv))   ystr_csv   <- "YSTR_table.csv"
if (is.null(haplogroup_hierarchy_csv)) haplogroup_hierarchy_csv <- "../input/Haplogroup_Hierarchy.csv"
if (is.null(root_out))   root_out <- "../output"

# --- Main switches & parameters ---
if (!exists("HH_range"))               HH_range     <- 1:10
if (!exists("use_equalize"))           use_equalize <- FALSE
if (!exists("thr_common"))             thr_common   <- 0.05
if (!exists("thr_rare"))               thr_rare     <- 0.01
if (!exists("jsd_eps"))                jsd_eps      <- 1e-12
if (!exists("make_multilayer_sankey")) make_multilayer_sankey <- TRUE
if (!exists("multilayer_max_k"))       multilayer_max_k <- 4
if (!exists("plot_base_height"))       plot_base_height <- 6
if (!exists("plot_base_width"))        plot_base_width  <- 8

# Multilayer alluvial display controls
if (!exists("perhap_max_plot")) perhap_max_plot <- 25
if (!exists("perhap_min_freq")) perhap_min_freq <- 0.005

# Sankey sizes
if (!exists("nature_sankey_width"))  nature_sankey_width  <- 8.5
if (!exists("nature_sankey_height")) nature_sankey_height <- 5.0

# Bootstrap controls
if (!exists("bootstrap_enable"))       bootstrap_enable       <- TRUE
if (!exists("bootstrap_n_per_sample")) bootstrap_n_per_sample <- 10000 # 10000
if (!exists("bootstrap_n_reps"))       bootstrap_n_reps       <- 1000 # 1000
if (!exists("bootstrap_pop_levels"))   bootstrap_pop_levels   <- NULL
if (!exists("bootstrap_ncol_facets"))  bootstrap_ncol_facets  <- 10
if (!exists("bootstrap_per_col_w"))    bootstrap_per_col_w    <- 2
if (!exists("bootstrap_per_row_h"))    bootstrap_per_row_h    <- 2
if (!exists("bootstrap_keep_only_levels")) bootstrap_keep_only_levels <- FALSE

# Significance annotation controls for bootstrap freq plots
if (!exists("freqplot_show_sig"))     freqplot_show_sig     <- FALSE  # set TRUE to show brackets
if (!exists("freqplot_sig_alpha"))    freqplot_sig_alpha    <- 0.05   # adjusted p-value threshold
if (!exists("freqplot_sig_ref_pop"))  freqplot_sig_ref_pop  <- NULL   # e.g., "Early_East_Asia"; if NULL, use top-N by p
if (!exists("freqplot_sig_max_pairs"))freqplot_sig_max_pairs<- 10     # max bracket count per hap

# Quasirandom jitter defaults
qr_method <- "pseudorandom"
qr_width  <- 0.22
qr_alpha  <- 0.60
qr_size   <- 0.9

dir.create(root_out, showWarnings = FALSE, recursive = TRUE)

# --- Dependencies ---
need_pkgs <- c(
  "readxl","dplyr","tidyr","stringr","ggplot2","ggrepel",
  "FactoMineR","vegan","pheatmap","ggalluvial","scales","ape","grid",
  "rstatix","readr","ggbeeswarm","ggsignif"
)
for (p in need_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    install.packages(p, repos = "https://cloud.r-project.org")
  }
}
library(readxl); library(dplyr); library(tidyr); library(stringr)
library(ggplot2); library(ggrepel); library(FactoMineR)
library(vegan); library(pheatmap); library(ggalluvial); library(scales); library(ape); library(grid)
library(rstatix); library(readr); library(ggbeeswarm); library(ggsignif)

suppressWarnings({
  ok <- try(Sys.setlocale("LC_CTYPE", "en_US.UTF-8"), silent = TRUE)
  if (inherits(ok, "try-error") || is.na(ok)) {
    if (.Platform$OS.type == "windows") {
      try(Sys.setlocale("LC_CTYPE","Chinese_China.936"), silent=TRUE)
    }
  }
})
options(encoding = "UTF-8")

`%||%` <- function(a, b) if (!is.null(a)) a else b
message2 <- function(fmt, ...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), sprintf(fmt, ...), "\n")
sanitize_label_ascii <- function(x){ x <- as.character(x); y <- iconv(x, from = "", to = "ASCII//TRANSLIT", sub = "_"); y[is.na(y)] <- "_"; y }
gsub_bytes <- function(pattern, replacement, x) gsub(pattern, replacement, x, useBytes = TRUE)

# --- Input helpers ---
safe_read_input <- function(xlsx){
  df <- readxl::read_xlsx(xlsx)
  cn <- tolower(colnames(df))
  id_idx  <- which(cn %in% c("sample_id","sample","id"))[1]
  hap_idx <- which(cn %in% c("haplogroup_full","haplogroup","hap"))[1]
  pop_idx <- which(cn %in% c("population","pop","group"))[1]
  if (any(is.na(c(id_idx, hap_idx, pop_idx)))) {
    if (length(cn) >= 3 && all(grepl("^v[0-9]+$", cn[1:3], ignore.case = TRUE))) {
      colnames(df)[1:3] <- c("sample_id","haplogroup_full","population")
    } else {
      if (ncol(df) < 3) stop("Input Excel must contain at least: sample_id, haplogroup_full, population")
      colnames(df)[1:3] <- c("sample_id","haplogroup_full","population")
    }
  } else { df <- df[, c(id_idx, hap_idx, pop_idx)]; colnames(df) <- c("sample_id","haplogroup_full","population") }
  df
}
sanitize_hap <- function(x){ x <- stringr::str_replace_all(x, "[^A-Za-z0-9]", ""); x <- tolower(x); if (nchar(x[1]) > 0) substr(x,1,1) <- toupper(substr(x,1,1)); x }
hap_substr <- function(h, HH) substr(h, 1, HH)

# --- NEW: read haplogroup hierarchy from file (wide format) ---
# Expected format: sample_id, hap_level_1, hap_level_2, ..., hap_level_N
read_haplogroup_hierarchy <- function(filepath){
  if (!file.exists(filepath)) return(NULL)
  df <- tryCatch(
    readr::read_csv(filepath, show_col_types = FALSE, col_types = cols(.default = col_character())),
    error = function(e) NULL
  )
  if (is.null(df) || !("sample_id" %in% colnames(df))) {
    message2("Haplogroup hierarchy file not found or missing 'sample_id' column: %s", filepath)
    return(NULL)
  }
  # Identify hap_level_* columns (or any numeric-suffix columns after sample_id)
  hap_cols <- grep("^hap_level_|^hap[0-9]", colnames(df), ignore.case = TRUE, value = TRUE)
  if (length(hap_cols) == 0) {
    message2("No hap_level_* columns found in %s; falling back to character truncation mode.", filepath)
    return(NULL)
  }
  df <- df[, c("sample_id", hap_cols), drop = FALSE]
  # Sort hap columns numerically if they follow pattern hap_level_1, hap_level_2, ...
  hap_nums <- as.numeric(gsub(".*_([0-9]+)$", "\\1", hap_cols))
  if (all(!is.na(hap_nums))) {
    idx <- order(hap_nums)
    hap_cols <- hap_cols[idx]
    df <- df[, c("sample_id", hap_cols), drop = FALSE]
  }
  colnames(df)[2:ncol(df)] <- paste0("hap_level_", seq_len(length(hap_cols)))
  message2("Loaded haplogroup hierarchy with %d levels from: %s", length(hap_cols), filepath)
  df
}

# --- NEW: merge hierarchy into main dataframe ---
merge_haplogroup_hierarchy <- function(df_main, df_hierarchy){
  if (is.null(df_hierarchy)) return(df_main)
  df_main <- dplyr::left_join(df_main, df_hierarchy, by = "sample_id")
  message2("Merged haplogroup hierarchy: %d samples matched", sum(!is.na(df_main$hap_level_1)))
  df_main
}

freq_table <- function(df, hap_col="hap", by_col="population"){
  df %>% count(.data[[by_col]], .data[[hap_col]], name="n") %>% group_by(.data[[by_col]]) %>% mutate(freq = n/sum(n)) %>% ungroup()
}
make_freq_matrix <- function(freq_df, by_col="population", hap_col="hap", val_col="freq"){
  freq_df %>% select(all_of(c(by_col,hap_col,val_col))) %>%
    tidyr::pivot_wider(names_from = all_of(hap_col), values_from = all_of(val_col), values_fill = 0) %>% as.data.frame()
}

# --- Diversity & distances ---
shannon <- function(p){ p <- p[p>0]; -sum(p*log(p)) }
simpson <- function(p){ 1 - sum(p^2) }
heff_shannon <- function(p){ exp(shannon(p)) }
heff_simpson <- function(p){ 1/sum(p^2) }
TVD <- function(p,q){ 0.5*sum(abs(p-q)) }
JSD_pair <- function(p,q, eps=jsd_eps){ p <- p + eps; p <- p/sum(p); q <- q + eps; q <- q/sum(q); m <- 0.5*(p+q); KL <- function(a,b) sum(a*log(a/b)); sqrt(0.5*KL(p,m) + 0.5*KL(q,m)) }
pairwise_matrix <- function(mat, FUN){ n <- nrow(mat); D <- matrix(0, n, n); for (i in seq_len(n)) { D[i,i] <- 0; if (i < n) { for (j in (i+1):n) { v <- FUN(mat[i,], mat[j,]); D[i,j] <- v; D[j,i] <- v }}}; rownames(D) <- rownames(mat); colnames(D) <- rownames(mat); D }
.nei_I <- function(p,q){ num <- sum(p*q); den <- sqrt(sum(p^2)*sum(q^2)); if (den==0) return(NA_real_); num/den }
.nei_D <- function(p,q){ I <- .nei_I(p,q); if (is.na(I) || I<=0) return(NA_real_); -log(I) }
.DA_pair <- function(p,q){ 1 - sum(sqrt(p*q)) }
pairwise_neiD <- function(mat){ pairwise_matrix(mat, .nei_D) }
pairwise_DA   <- function(mat){ pairwise_matrix(mat, .DA_pair) }

# --- Y-STR Rst ---
compute_pairwise_Rst <- function(ystr_df, pop_col="population"){
  stopifnot(pop_col %in% names(ystr_df))
  loci <- setdiff(colnames(ystr_df), c("sample_id", pop_col))
  loci <- loci[sapply(ystr_df[loci], function(x) is.numeric(x) || is.integer(x))]
  pops <- unique(ystr_df[[pop_col]]); P <- length(pops)
  R <- matrix(NA_real_, P, P, dimnames=list(pops,pops))
  rst_two <- function(xA,xB){
    xA <- xA[is.finite(xA)]; xB <- xB[is.finite(xB)]
    if (length(xA)<2 || length(xB)<2) return(NA_real_)
    nA <- length(xA); nB <- length(xB); nT <- nA+nB
    mA <- mean(xA); mB <- mean(xB); mT <- (nA*mA + nB*mB)/nT
    SSB <- nA*(mA-mT)^2 + nB*(mB-mT)^2
    SSW <- sum((xA-mA)^2) + sum((xB-mB)^2)
    if ((SSB+SSW)<=0) return(NA_real_)
    SSB/(SSB+SSW)
  }
  for (i in seq_len(P)){
    for (j in seq_len(P)){
      if (j<i) { R[i,j] <- R[j,i]; next }
      if (i==j) { R[i,j] <- 0; next }
      rloc <- c()
      for (L in loci){
        xi <- ystr_df[ystr_df[[pop_col]]==pops[i], L][[1]]
        xj <- ystr_df[ystr_df[[pop_col]]==pops[j], L][[1]]
        rloc <- c(rloc, rst_two(xi,xj))
      }
      if (length(rloc)) R[i,j] <- mean(rloc, na.rm = TRUE)
      if (is.nan(R[i,j])) R[i,j] <- NA_real_
    }
  }
  R
}

# --- Enrichment (Fisher) ---
hap_enrichment <- function(df, hap_col="hap", pop_col="population"){
  all_pops <- unique(df[[pop_col]]); all_haps <- unique(df[[hap_col]]); res <- list()
  for (h in all_haps){
    df$flag <- df[[hap_col]]==h
    for (g in all_pops){
      a <- sum(df$flag & df[[pop_col]]==g); b <- sum(!df$flag & df[[pop_col]]==g)
      c <- sum(df$flag & df[[pop_col]]!=g); d <- sum(!df$flag & df[[pop_col]]!=g)
      mat <- matrix(c(a,b,c,d), nrow=2, byrow=TRUE)
      if (all(mat>=0)){ ft <- suppressWarnings(fisher.test(mat))
        res[[length(res)+1]] <- data.frame(hap=h, population=g, a=a,b=b,c=d[1],d=d[2], p=ft$p.value, odds=unname(ft$estimate)) }
    }
  }
  out <- do.call(rbind, res); out$padj <- p.adjust(out$p, method="BH"); out[order(out$padj, out$p),]
}

# --- Meta helper ---
safe_group_from_meta <- function(meta_in){
  if (is.null(meta_in) || !is.data.frame(meta_in) || ncol(meta_in)==0) return(NULL)
  meta_lc <- meta_in; colnames(meta_lc) <- tolower(colnames(meta_lc))
  if (!("population" %in% colnames(meta_lc))){
    cn <- tolower(names(meta_lc)); pid <- which(cn %in% c("population","pop","group_name"))
    if (length(pid)==1) names(meta_lc)[pid] <- "population"
  }
  if ("group" %in% colnames(meta_lc)){ out <- meta_lc[, c("population","group"), drop=FALSE] }
  else if ("region" %in% colnames(meta_lc)){ out <- meta_lc[, c("population","region"), drop=FALSE]; names(out)[2] <- "group" }
  else return(NULL)
  out$population <- sanitize_label_ascii(out$population); out$group <- sanitize_label_ascii(out$group); out
}

# --- Distance matrix sanitization ---
sanitize_dist_matrix <- function(D, name = "D", eps = 1e-6) {
  if (is.null(D)) stop("distance matrix is NULL: ", name)
  D <- as.matrix(D)
  if (!isSymmetric(D, tol = 1e-8)) D <- (D + t(D)) / 2
  if (any(!is.finite(D))) {
    finite_vals <- D[is.finite(D)]; repl <- if (length(finite_vals)) max(finite_vals) + eps else 1
    D[!is.finite(D)] <- repl
  }
  diag(D) <- 0; up <- D[upper.tri(D)]
  if (length(up) > 1) {
    if (all(abs(up - up[1]) < 1e-12)) {
      n <- nrow(D); if (n >= 2) { for (i in seq_len(n)) for (j in seq_len(n)) if (i != j) D[i,j] <- abs(i-j)*eps; diag(D) <- 0 }
    }
  }
  return(D)
}

# --- Colors & ordering ---
.okabe_ito <- c("#0072B2","#56B4E9","#009E73","#E69F00","#CC79A7","#D55E00","#999999","#F0E442")
.nature_seeds <- c("#4477AA","#66CCEE","#228833","#CCBB44","#EE6677","#AA3377",
                   "#332288","#88CCEE","#44AA99","#117733","#999933","#DDCC77",
                   "#CC6677","#882255","#AA4499","#DDDDDD",
                   "#1F77B4","#FF7F0E","#2CA02C","#D62728","#9467BD","#8C564B",
                   "#E377C2","#7F7F7F","#BCBD22","#17BECF")
.seeds <- unique(c(.okabe_ito, .nature_seeds))
.hex_to_lab <- function(hex){ hex <- as.character(hex); rgb <- t(grDevices::col2rgb(hex)/255); lab <- grDevices::convertColor(rgb, from="sRGB", to="Lab", scale.in=1); lab <- as.matrix(lab); colnames(lab) <- c("L","a","b"); lab }
.gen_candidates <- function(){ levL <- c(80,72,64,56,48); levC <- c(60,66,62,58,54); out <- character(); for (i in seq_along(levL)){ cols <- grDevices::hcl(h = seq(0, 359, by = 1.5), c = levC[i], l = levL[i]); out <- c(out, cols) }; unique(out) }
glasbey_like <- function(n, seed_colors = .seeds, max_n = 2048){
  n <- min(n, max_n); seed_colors <- unique(seed_colors)
  if (n <= length(seed_colors)) return(seed_colors[seq_len(n)])
  cand <- setdiff(.gen_candidates(), seed_colors); sel  <- seed_colors
  sel_lab  <- .hex_to_lab(sel); cand_lab <- .hex_to_lab(cand)
  set.seed(42); ord <- sample.int(nrow(cand_lab)); cand_lab <- cand_lab[ord, , drop=FALSE]; cand <- cand[ord]
  min_dist_to_sel <- function(row, sel_lab){ min(sqrt(rowSums((t(sel_lab) - row)^2))) }
  dist_min <- apply(cand_lab, 1, min_dist_to_sel, sel_lab = sel_lab)
  while (length(sel) < n && nrow(cand_lab) > 0){
    pick <- which.max(dist_min); sel <- c(sel, cand[pick]); new_row <- cand_lab[pick, , drop=FALSE]
    sel_lab <- rbind(sel_lab, new_row)
    dnew <- sqrt(rowSums((cand_lab - matrix(new_row, nrow=nrow(cand_lab), ncol=3, byrow=TRUE))^2))
    dist_min <- pmin(dist_min, dnew); cand_lab <- cand_lab[-pick, , drop=FALSE]; cand <- cand[-pick]; dist_min <- dist_min[-pick]
  }
  sel[seq_len(n)]
}
get_palette_for_levels <- function(fac_levels, min_n = 0){
  lv <- as.character(unique(fac_levels)); n  <- max(length(lv), min_n)
  pal <- glasbey_like(n); setNames(pal[seq_len(length(lv))], lv)
}
nature_hm_seq <- function(n = 100){
  grDevices::colorRampPalette(c("#F7FBFF","#DEEBF7","#C6DBEF",
                                "#9ECAE1","#6BAED6","#4292C6",
                                "#2171B5","#08519C","#08306B"))(n)
}
save_png_pdf <- function(p, file_path_noext, width, height){
  suppressWarnings(ggsave(paste0(file_path_noext, ".png"), p, width=width, height=height, dpi=400, limitsize = FALSE))
  suppressWarnings(ggsave(paste0(file_path_noext, ".pdf"), p, width=width, height=height, limitsize = FALSE))
}
order_pop_by_dominant_hap <- function(freq_mat_df){
  long <- freq_mat_df %>% rename(population = 1) %>% pivot_longer(-population, names_to = "hap", values_to = "freq")
  long$hap <- droplevels(factor(long$hap))
  pop_top <- long %>% group_by(population) %>% slice_max(freq, n = 1, with_ties = FALSE) %>% ungroup() %>% rename(top_hap = hap, top_freq = freq)
  hap_rank <- pop_top %>% group_by(top_hap) %>% summarise(max_freq_any_pop = max(top_freq, na.rm = TRUE), .groups = "drop") %>% arrange(desc(max_freq_any_pop)) %>% mutate(hap_order = row_number())
  pop_order_df <- pop_top %>% left_join(hap_rank, by = "top_hap") %>% arrange(hap_order, desc(top_freq), population)
  pop_order_df$population
}

# --- Plotting ---
plot_distance_heatmap <- function(D, out_file, title="Distance heatmap"){
  D <- sanitize_dist_matrix(D, name = title); cols <- nature_hm_seq(100); base <- tools::file_path_sans_ext(out_file)
  hc <- hclust(as.dist(D), method = "ward.D2"); ord <- hc$order; D2  <- D[ord, ord, drop = FALSE]
  tryCatch({
    png(paste0(base, ".png"), width  = plot_base_width*100, height = plot_base_height*100)
    suppressWarnings(pheatmap::pheatmap(D2, color=cols, cluster_rows=FALSE, cluster_cols=FALSE, main=title,
                                        border_color=NA, legend=TRUE, show_rownames=TRUE, show_colnames=TRUE))
    dev.off()
  }, error = function(e) { tryCatch(dev.off(), error = function(e2) {}) })
  tryCatch({
    pdf(paste0(base, ".pdf"), width = plot_base_width, height = plot_base_height)
    suppressWarnings(pheatmap::pheatmap(D2, color=cols, cluster_rows=FALSE, cluster_cols=FALSE, main=title,
                                        border_color=NA, legend=TRUE, show_rownames=TRUE, show_colnames=TRUE))
    dev.off()
  }, error = function(e) { tryCatch(dev.off(), error = function(e2) {}) })
}
plot_dendrogram <- function(D, out_file){
  D <- sanitize_dist_matrix(D, name = out_file); hc <- hclust(as.dist(D), method="ward.D2"); base <- tools::file_path_sans_ext(out_file)
  tryCatch({
    pdf(paste0(base, ".pdf"), width=plot_base_width, height=plot_base_height)
    plot(hc, xlab="", sub="", main="Hierarchical clustering (ward.D2)")
    dev.off()
  }, error = function(e) { tryCatch(dev.off(), error = function(e2) {}) })
  tryCatch({
    png(paste0(base, ".png"), width=plot_base_width*100, height=plot_base_height*100)
    plot(hc, xlab="", sub="", main="Hierarchical clustering (ward.D2)")
    dev.off()
  }, error = function(e) { tryCatch(dev.off(), error = function(e2) {}) })
  phy <- ape::as.phylo(hc); phy$tip.label <- sanitize_label_ascii(gsub_bytes("[^A-Za-z0-9_.-]", "_", phy$tip.label))
  ape::write.tree(phy, file=paste0(base, ".nwk")); message2("Exported Newick tree: %s", normalizePath(paste0(base, ".nwk")))
}
plot_nj_tree <- function(D, out_file, title="Neighbor-Joining tree"){
  D <- sanitize_dist_matrix(D, name = out_file); tr <- ape::nj(as.dist(D)); base <- tools::file_path_sans_ext(out_file)
  tryCatch({
    pdf(paste0(base, ".pdf"), width=plot_base_width, height=plot_base_height)
    plot(tr, main=title, cex=0.8)
    dev.off()
  }, error = function(e) { tryCatch(dev.off(), error = function(e2) {}) })
  tryCatch({
    png(paste0(base, ".png"), width=plot_base_width*100, height=plot_base_height*100)
    plot(tr, main=title, cex=0.8)
    dev.off()
  }, error = function(e) { tryCatch(dev.off(), error = function(e2) {}) })
  tr$tip.label <- sanitize_label_ascii(gsub_bytes("[^A-Za-z0-9_.-]", "_", tr$tip.label))
  ape::write.tree(tr, file=paste0(base, ".nwk")); message2("Exported NJ Newick: %s", normalizePath(paste0(base, ".nwk")))
}
plot_stacked_bars <- function(freq_mat_df, out_file){
  long <- freq_mat_df %>% rename(population = 1) %>% pivot_longer(-population, names_to="hap", values_to="freq")
  long$hap <- droplevels(factor(long$hap))
  pal <- get_palette_for_levels(levels(long$hap), min_n = max(256, nlevels(long$hap)))
  p <- ggplot(long, aes(x=population, y=freq, fill=hap)) +
    geom_bar(stat="identity", width=0.8) +
    scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
    scale_fill_manual(values=pal, limits=names(pal), drop=FALSE, na.value="#BDBDBD") +
    theme_bw(base_size=12) +
    theme(axis.text.x = element_text(angle=60, hjust=1),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank()) +
    labs(x="Population", y="Frequency", fill="Haplogroup")
  base <- tools::file_path_sans_ext(out_file)
  save_png_pdf(p, base, width=plot_base_width*1.6, height=plot_base_height*1.2)
}
run_pca_plot <- function(freq_mat_df, out_png, out_scores_csv){
  mat <- as.matrix(freq_mat_df[,-1]); rownames(mat) <- sanitize_label_ascii(freq_mat_df[[1]])
  pca <- prcomp(mat, center=TRUE, scale.=FALSE); var_exp <- (pca$sdev^2)/sum(pca$sdev^2)
  scores <- data.frame(population=rownames(pca$x), pca$x[,1:5, drop=FALSE], check.names=FALSE)
  grp <- safe_group_from_meta(if (exists("meta_in", inherits=TRUE)) meta_in else NULL)
  if (!is.null(grp)) { scores <- dplyr::left_join(scores, grp, by="population") } else {
    hc <- hclust(dist(mat), method="ward.D2"); cl <- cutree(hc, k=min(4, nrow(mat))); scores$group <- factor(cl)
  }
  write.csv(scores, out_scores_csv, row.names=FALSE, fileEncoding="UTF-8")
  pc1 <- sprintf("PC1 (%.1f%%)", 100*var_exp[1]); pc2 <- sprintf("PC2 (%.1f%%)", 100*var_exp[2])
  cols <- get_palette_for_levels(levels(factor(scores$group)))
  p <- ggplot(scores, aes(PC1, PC2, label=population, color=group)) +
    geom_point(size=3) + geom_text_repel(hjust=0, max.overlaps=100) +
    scale_color_manual(values=cols, limits=names(cols), drop=FALSE, na.value="#BDBDBD") +
    theme_bw(base_size=12) + theme(panel.grid = element_blank()) +
    labs(title="PCA of haplogroup frequencies", x=pc1, y=pc2, color="Group")
  base <- tools::file_path_sans_ext(out_png); save_png_pdf(p, base, width=plot_base_width, height=plot_base_height)
}
run_ca_plot <- function(count_mat_df, out_png, out_coords_csv){
  mat <- as.data.frame(count_mat_df); rownames(mat) <- sanitize_label_ascii(mat[[1]]); mat <- mat[,-1, drop=FALSE]
  ca <- FactoMineR::CA(mat, graph=FALSE)
  row_coords <- data.frame(population=rownames(ca$row$coord), ca$row$coord, check.names=FALSE)
  col_coords <- data.frame(hap=rownames(ca$col$coord), ca$col$coord, check.names=FALSE)
  eig <- as.data.frame(ca$eig)
  row_contrib <- data.frame(population=rownames(ca$row$contrib), ca$row$contrib, check.names=FALSE)
  col_contrib <- data.frame(hap=rownames(ca$col$contrib), ca$col$contrib, check.names=FALSE)
  write.csv(row_coords, sub("coords","row_coords",out_coords_csv), row.names=FALSE, fileEncoding="UTF-8")
  write.csv(col_coords, sub("coords","col_coords",out_coords_csv), row.names=FALSE, fileEncoding="UTF-8")
  write.csv(eig,        sub("coords","eigenvalues",out_coords_csv), row.names=FALSE, fileEncoding="UTF-8")
  write.csv(row_contrib,sub("coords","row_contrib",out_coords_csv), row.names=FALSE, fileEncoding="UTF-8")
  write.csv(col_contrib,sub("coords","col_contrib",out_coords_csv), row.names=FALSE, fileEncoding="UTF-8")
  scores <- row_coords; var_exp <- eig[,2]
  grp <- safe_group_from_meta(if (exists("meta_in", inherits=TRUE)) meta_in else NULL)
  if (!is.null(grp)) { scores <- dplyr::left_join(scores, grp, by="population") } else {
    hc <- hclust(dist(as.matrix(mat)), method="ward.D2"); cl <- cutree(hc, k=min(4, nrow(mat))); scores$group <- factor(cl)
  }
  rows12 <- scores[, c("population","Dim 1","Dim 2","group")]; colnames(rows12) <- c("population","Dim1","Dim2","group")
  cols12 <- col_coords[, c("hap","Dim 1","Dim 2")]; colnames(cols12) <- c("hap","Dim1","Dim2")
  pc1 <- sprintf("Dim1 (%.1f%%)", var_exp[1]); pc2 <- sprintf("Dim2 (%.1f%%)", var_exp[2])
  cols_grp <- get_palette_for_levels(levels(factor(rows12$group)))
  p <- ggplot() +
    geom_point(data=rows12, aes(Dim1, Dim2, color=group), size=3) +
    geom_text_repel(data=rows12, aes(Dim1, Dim2, label=population, color=group),
                    hjust=0, max.overlaps=100) +
    geom_segment(data=cols12, aes(x=0,y=0,xend=Dim1,yend=Dim2),
                 arrow=arrow(length=unit(0.02,"npc")), alpha=0.6) +
    geom_text_repel(data=cols12, aes(Dim1, Dim2, label=hap), size=3) +
    scale_color_manual(values=cols_grp, limits=names(cols_grp), drop=FALSE, na.value="#BDBDBD") +
    theme_bw(base_size=12) + theme(panel.grid = element_blank()) +
    labs(title="CA biplot (rows=populations, cols=haplogroups)", x=pc1, y=pc2, color="Group")
  base <- tools::file_path_sans_ext(out_png); save_png_pdf(p, base, width=plot_base_width, height=plot_base_height)
}
run_mds_plot <- function(D, out_png){
  D <- sanitize_dist_matrix(D, name = out_png)
  conf <- cmdscale(as.dist(D), k=2, eig=TRUE)
  pts <- as.data.frame(conf$points); pts$population <- sanitize_label_ascii(rownames(pts))
  colnames(pts)[1:2] <- c("Dim1","Dim2")
  eig <- conf$eig; pos <- eig[eig>0]; var_exp <- if (length(pos)) pos/sum(pos) else c(NA_real_, NA_real_)
  var_lab <- c(sprintf("Dim1 (%.1f%%)", 100*ifelse(length(var_exp)>=1, var_exp[1], NA)),
               sprintf("Dim2 (%.1f%%)", 100*ifelse(length(var_exp)>=2, var_exp[2], NA)))
  grp <- safe_group_from_meta(if (exists("meta_in", inherits=TRUE)) meta_in else NULL)
  if (!is.null(grp)) { pts <- dplyr::left_join(pts, grp, by="population") } else {
    hc <- hclust(as.dist(D), method="ward.D2"); cl <- cutree(hc, k=min(4,nrow(D))); pts$group <- factor(cl)
  }
  base_csv <- tools::file_path_sans_ext(out_png)
  write.csv(pts, paste0(base_csv, "_MDS_points.csv"), row.names=FALSE, fileEncoding="UTF-8")
  write.csv(data.frame(eigenvalue=eig, component=seq_along(eig)),
            paste0(base_csv, "_MDS_eigenvalues.csv"), row.names=FALSE, fileEncoding="UTF-8")
  cols <- get_palette_for_levels(levels(factor(pts$group)))
  p <- ggplot(pts, aes(Dim1, Dim2, label=population, color=group)) +
    geom_point(size=3) + geom_text_repel(hjust=0, max.overlaps=100) +
    scale_color_manual(values=cols, limits=names(cols), drop=FALSE, na.value="#BDBDBD") +
    theme_bw(base_size=12) + theme(panel.grid = element_blank()) +
    labs(title="MDS (cmdscale) on distance matrix", x=var_lab[1], y=var_lab[2], color="Group")
  base <- tools::file_path_sans_ext(out_png); save_png_pdf(p, base, width=plot_base_width, height=plot_base_height)
}

# --- Simple alluvials ---
plot_alluvial_simple <- function(df, hap_col, pop_col, out_png){
  tab <- df %>% count(.data[[hap_col]], .data[[pop_col]], name="n") %>% rename(hap=1, population=2)
  tab$hap <- droplevels(factor(tab$hap)); tab$population <- factor(sanitize_label_ascii(tab$population))
  pal <- get_palette_for_levels(levels(tab$hap), min_n = max(512, nlevels(tab$hap)))
  p <- ggplot(tab, aes(y=n, axis1=hap, axis2=population)) +
    geom_alluvium(aes(fill=hap), width=1/12, knot.pos=0.4, alpha=0.9) +
    geom_stratum(aes(fill=hap), width=1/12, color="grey30") +
    geom_text(stat="stratum", aes(label=after_stat(stratum)), hjust=0.5, vjust=0.5, size=3) +
    scale_x_discrete(limits=c("Haplogroup","Population"), expand=c(.1,.05)) +
    scale_fill_manual(values=pal, limits=names(pal), drop=TRUE, na.value="#BDBDBD") +
    theme_bw(base_size=12) + theme(panel.grid = element_blank()) +
    labs(y="Count", x="", fill="Haplogroup", title="Alluvial: Haplogroup → Population (Counts)")
  base <- tools::file_path_sans_ext(out_png); save_png_pdf(p, base, width=plot_base_width*1.8, height=plot_base_height*1.3)
}
plot_alluvial_simple_freq <- function(freq_mat_df, out_png){
  long <- freq_mat_df %>% rename(population=1) %>% pivot_longer(-population, names_to="hap", values_to="freq")
  long$hap <- droplevels(factor(long$hap)); long$population <- factor(sanitize_label_ascii(long$population))
  pal <- get_palette_for_levels(levels(long$hap), min_n = max(512, nlevels(long$hap)))
  p <- ggplot(long, aes(y=freq, axis1=hap, axis2=population)) +
    geom_alluvium(aes(fill=hap), width=1/12, knot.pos=0.4, alpha=0.9) +
    geom_stratum(aes(fill=hap), width=1/12, color="grey30") +
    geom_text(stat="stratum", aes(label=after_stat(stratum)), hjust=0.5, vjust=0.5, size=3) +
    scale_x_discrete(limits=c("Haplogroup","Population"), expand=c(.1,.05)) +
    scale_y_continuous(labels=scales::percent_format(accuracy=1)) +
    scale_fill_manual(values=pal, limits=names(pal), drop=TRUE, na.value="#BDBDBD") +
    theme_bw(base_size=12) + theme(panel.grid = element_blank()) +
    labs(y="Frequency", x="", fill="Haplogroup", title="Alluvial: Haplogroup → Population (Frequency)")
  base <- tools::file_path_sans_ext(out_png); save_png_pdf(p, base, width=plot_base_width*1.8, height=plot_base_height*1.3)
}

# --- Multilayer alluvial core (patched) ---
make_multilayer_long <- function(df, HH_chain = c(1,2,3,4)){
  tmp <- df
  # NEW: Use hap_level_* columns directly instead of truncating hap_clean
  for (k in HH_chain) {
    hap_col <- paste0("hap_level_", k)
    if (hap_col %in% colnames(tmp)) {
      tmp[[paste0("hap", k)]] <- tmp[[hap_col]]
    } else {
      warning(sprintf("Column %s not found; skipping level %d in multilayer", hap_col, k))
    }
  }
  axes <- c(paste0("hap", HH_chain), "population")
  # Filter axes to only those that exist in tmp
  axes_exist <- axes[axes %in% c(colnames(tmp), "population")]
  if (length(axes_exist) == 0) stop("No valid axes found for multilayer alluvial")
  
  long <- tmp %>% dplyr::count(dplyr::across(dplyr::all_of(axes_exist)), name = "n") %>%
    dplyr::group_by(population) %>% dplyr::mutate(freq = n / sum(n)) %>% dplyr::ungroup()
  for (ax in axes_exist) long[[ax]] <- as.character(long[[ax]])
  stopifnot(all(vapply(long[axes_exist], is.character, TRUE))); long$population <- as.character(long$population)
  deepest_col    <- axes_exist[length(axes_exist) - 1]; deepest_levels <- sort(unique(long[[deepest_col]]))
  pal_deep <- get_palette_for_levels(deepest_levels, min_n = max(1024, length(deepest_levels)))
  lighten_col <- function(col, f = 0.3){ v <- grDevices::col2rgb(col)/255; v <- v + (1 - v) * f; rgb(v[1], v[2], v[3]) }
  hap_axes <- axes_exist[-length(axes_exist)]; n_hap_axes <- length(hap_axes)
  pal_all <- pal_deep
  if (n_hap_axes > 1) {
    for (i in seq_len(n_hap_axes - 1)) {
      ax <- hap_axes[i]; levs_ax <- sort(unique(long[[ax]]))
      for (lv in levs_ax) {
        cand <- deepest_levels[startsWith(deepest_levels, lv)]
        base_col <- if (!length(cand)) "#BDBDBD" else if (length(cand) == 1) { pal_deep[cand] } else {
          cols <- pal_deep[cand]; rgbs <- sapply(cols, function(z) grDevices::col2rgb(z)/255); m <- rowMeans(rgbs); rgb(m[1], m[2], m[3])
        }
        frac <- (n_hap_axes - i)/(n_hap_axes - 1 + 1e-9); pal_all[[lv]] <- lighten_col(base_col, f = 0.35 * frac)
      }
    }
  }
  pop_levels <- sort(unique(long$population)); pal_all <- c(pal_all, setNames(rep("#E0E0E0", length(pop_levels)), pop_levels))
  list(long=long, axes=axes_exist, deepest=deepest_col, deepest_levels=deepest_levels, pal_all=pal_all, pop_levels=pop_levels)
}
.compute_label_keep <- function(long, hap_axes, label_min_freq = 0.01){
  keep <- character(0)
  for (ax in hap_axes){
    tmp <- long %>% dplyr::group_by(population, lvl = .data[[ax]]) %>%
      dplyr::summarise(f = sum(freq, na.rm = TRUE), .groups = "drop") %>%
      dplyr::group_by(lvl) %>% dplyr::summarise(avg_freq = mean(f, na.rm = TRUE), .groups = "drop")
    keep <- union(keep, as.character(tmp$lvl[tmp$avg_freq >= label_min_freq]))
  }
  keep[!is.na(keep) & nzchar(keep)]
}
plot_alluvial_multilayer_counts_from_long <- function(mdat, out_png, label_min_freq = 0.01){
  long <- mdat$long; axes <- mdat$axes; deepest <- mdat$deepest; pal_all <- mdat$pal_all
  hap_axes <- axes[-length(axes)]; pop_levels <- mdat$pop_levels; deep_levels <- mdat$deepest_levels
  label_keep <- .compute_label_keep(long, hap_axes, label_min_freq)
  mapping <- do.call(ggplot2::aes, c(list(y = quote(n)), setNames(lapply(axes, as.name), paste0("axis", seq_along(axes)))))
  p <- ggplot2::ggplot(long, mapping) +
    ggalluvial::geom_alluvium(ggplot2::aes_string(fill = deepest), width = 1/14, knot.pos = 0.38, alpha = 0.85) +
    ggalluvial::geom_stratum(ggplot2::aes(fill = after_stat(stratum)), width = 1/14, color = "grey30", alpha = 0.98, show.legend = FALSE) +
    ggplot2::geom_text(stat = "stratum",
      ggplot2::aes(label = dplyr::if_else(after_stat(stratum) %in% !!pop_levels,
                after_stat(stratum),
                dplyr::if_else(after_stat(stratum) %in% !!label_keep, after_stat(stratum), NA_character_))),
      size = 3, vjust = 0.5, hjust = 0.5) +
    ggplot2::scale_x_discrete(limits = c(paste0("HH=", sub("hap","", hap_axes)), "Population"), expand = c(.1,.05)) +
    ggplot2::scale_fill_manual(values = pal_all, breaks = deep_levels, drop = TRUE, na.value = "#BDBDBD") +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Haplogroup")) +
    ggplot2::theme_bw(base_size = 12) + ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 10), plot.margin = ggplot2::margin(10,40,10,10)) +
    ggplot2::coord_cartesian(clip = "off") + ggplot2::labs(y = "Count", x = "", title = "Multilayer Alluvial (Counts)")
  base <- tools::file_path_sans_ext(out_png)
  suppressWarnings(ggplot2::ggsave(paste0(base, ".png"), p, width = plot_base_width*2.0, height = plot_base_height*1.25, dpi = 400, limitsize=FALSE))
  suppressWarnings(ggplot2::ggsave(paste0(base, ".pdf"), p, width = plot_base_width*2.0, height = plot_base_height*1.25, limitsize=FALSE))
}
plot_alluvial_multilayer_freq_from_long <- function(mdat, out_png, label_min_freq = 0.01){
  long <- mdat$long; axes <- mdat$axes; deepest <- mdat$deepest; pal_all <- mdat$pal_all
  hap_axes <- axes[-length(axes)]; pop_levels <- mdat$pop_levels; deep_levels <- mdat$deepest_levels
  label_keep <- .compute_label_keep(long, hap_axes, label_min_freq)
  mapping <- do.call(ggplot2::aes, c(list(y = quote(freq)), setNames(lapply(axes, as.name), paste0("axis", seq_along(axes)))))
  p <- ggplot2::ggplot(long, mapping) +
    ggalluvial::geom_alluvium(ggplot2::aes_string(fill = deepest), width = 1/14, knot.pos = 0.38, alpha = 0.85) +
    ggalluvial::geom_stratum(ggplot2::aes(fill = after_stat(stratum)), width = 1/14, color = "grey30", alpha = 0.98, show.legend = FALSE) +
    ggplot2::geom_text(stat = "stratum",
      ggplot2::aes(label = dplyr::if_else(after_stat(stratum) %in% !!pop_levels,
                after_stat(stratum),
                dplyr::if_else(after_stat(stratum) %in% !!label_keep, after_stat(stratum), NA_character_))),
      size = 3, vjust = 0.5, hjust = 0.5) +
    ggplot2::scale_x_discrete(limits = c(paste0("HH=", sub("hap","", hap_axes)), "Population"), expand = c(.1,.05)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_fill_manual(values = pal_all, breaks = deep_levels, drop = TRUE, na.value = "#BDBDBD") +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Haplogroup")) +
    ggplot2::theme_bw(base_size = 12) + ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 10), plot.margin = ggplot2::margin(10,40,10,10)) +
    ggplot2::coord_cartesian(clip = "off") + ggplot2::labs(y = "Frequency", x = "", title = "Multilayer Alluvial (Frequency)")
  base <- tools::file_path_sans_ext(out_png)
  suppressWarnings(ggplot2::ggsave(paste0(base, ".png"), p, width = plot_base_width*2.0, height = plot_base_height*1.25, dpi = 400, limitsize=FALSE))
  suppressWarnings(ggplot2::ggsave(paste0(base, ".pdf"), p, width = plot_base_width*2.0, height = plot_base_height*1.25, limitsize=FALSE))
}
plot_alluvial_multilayer_freq_singlehap_raw <- function(mdat, root_hap, out_png, label_min_freq = 0.01){
  long <- mdat$long; axes <- mdat$axes; deepest <- mdat$deepest; pal_all <- mdat$pal_all; pop_levels <- mdat$pop_levels
  root_col <- axes[1]; sub_long <- long %>% dplyr::filter(.data[[root_col]] == root_hap)
  if (nrow(sub_long) == 0) { message(sprintf("Root %s has no data under this HH chain; skip raw plot", root_hap)); return(invisible(NULL)) }
  hap_axes <- axes[-length(axes)]; label_keep <- .compute_label_keep(sub_long, hap_axes, label_min_freq)
  mapping <- do.call(ggplot2::aes, c(list(y = quote(freq)), setNames(lapply(axes, as.name), paste0("axis", seq_along(axes)))))
  p <- ggplot2::ggplot(sub_long, mapping) +
    ggalluvial::geom_alluvium(ggplot2::aes_string(fill = deepest), width = 1/14, knot.pos = 0.38, alpha = 0.9) +
    ggalluvial::geom_stratum(ggplot2::aes(fill = as.character(after_stat(stratum))), width = 1/14, color = "grey30", show.legend = FALSE) +
    ggplot2::geom_text(stat  = "stratum",
      ggplot2::aes(label = ifelse(as.character(after_stat(stratum)) %in% !!pop_levels, as.character(after_stat(stratum)),
                                  ifelse(as.character(after_stat(stratum)) %in% !!label_keep, as.character(after_stat(stratum)), NA_character_)), hjust  = 0.5),
      size = 3, vjust = 0.5) +
    ggplot2::scale_x_discrete(limits = c(paste0("HH=", sub("hap","", hap_axes)), "Population"), expand = c(.1, .05)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_fill_manual(values = pal_all, breaks = unique(sub_long[[deepest]]), drop = TRUE, na.value = "#BDBDBD") +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Haplogroup")) +
    ggplot2::theme_bw(base_size = 12) + ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 10), plot.margin = ggplot2::margin(10,40,10,10)) +
    ggplot2::coord_cartesian(clip = "off") + ggplot2::labs(y = "Frequency", x = "", title = paste0("Multilayer Alluvial (", root_hap, ", raw)"))
  base <- tools::file_path_sans_ext(out_png)
  suppressWarnings(ggplot2::ggsave(paste0(base, ".png"), p, width = plot_base_width*2.0, height = plot_base_height*1.25, dpi = 400, limitsize=FALSE))
  suppressWarnings(ggplot2::ggsave(paste0(base, ".pdf"), p, width = plot_base_width*2.0, height = plot_base_height*1.25, limitsize=FALSE))
}
plot_alluvial_multilayer_freq_singlehap_scaled <- function(mdat, root_hap, out_png, label_min_freq = 0.01){
  long <- mdat$long; axes <- mdat$axes; deepest <- mdat$deepest; pal_all <- mdat$pal_all; pop_levels <- mdat$pop_levels
  root_col <- axes[1]; sub_long <- long %>% dplyr::filter(.data[[root_col]] == root_hap)
  if (nrow(sub_long) == 0) { message(sprintf("Root %s has no data under this HH chain; skip scaled plot", root_hap)); return(invisible(NULL)) }
  sub_long <- sub_long %>% dplyr::group_by(dplyr::across(dplyr::all_of(axes))) %>% dplyr::mutate(freq_scaled = freq / sum(freq)) %>% dplyr::ungroup()
  hap_axes <- axes[-length(axes)]; sub_long_lbl <- sub_long %>% dplyr::mutate(freq = freq_scaled); label_keep <- .compute_label_keep(sub_long_lbl, hap_axes, label_min_freq)
  mapping <- do.call(ggplot2::aes, c(list(y = quote(freq_scaled)), setNames(lapply(axes, as.name), paste0("axis", seq_along(axes)))))
  p <- ggplot2::ggplot(sub_long, mapping) +
    ggalluvial::geom_alluvium(ggplot2::aes_string(fill = deepest), width = 1/14, knot.pos = 0.38, alpha = 0.9) +
    ggalluvial::geom_stratum(ggplot2::aes(fill = as.character(after_stat(stratum))), width = 1/14, color = "grey30", show.legend = FALSE) +
    ggplot2::geom_text(stat  = "stratum",
      ggplot2::aes(label = ifelse(as.character(after_stat(stratum)) %in% !!pop_levels, as.character(after_stat(stratum)),
                                  ifelse(as.character(after_stat(stratum)) %in% !!label_keep, as.character(after_stat(stratum)), NA_character_)), hjust  = 0.5),
      size = 3, vjust = 0.5) +
    ggplot2::scale_x_discrete(limits = c(paste0("HH=", sub("hap","", hap_axes)), "Population"), expand = c(.1, .05)) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    ggplot2::scale_fill_manual(values = pal_all, breaks = unique(sub_long[[deepest]]), drop = TRUE, na.value = "#BDBDBD") +
    ggplot2::guides(fill = ggplot2::guide_legend(title = "Haplogroup")) +
    ggplot2::theme_bw(base_size = 12) + ggplot2::theme(panel.grid = ggplot2::element_blank(), axis.text.x = ggplot2::element_text(size = 10), plot.margin = ggplot2::margin(10,40,10,10)) +
    ggplot2::coord_cartesian(clip = "off") + ggplot2::labs(y = "Scaled frequency", x = "", title = paste0("Multilayer Alluvial (", root_hap, ", scaled)"))
  base <- tools::file_path_sans_ext(out_png)
  suppressWarnings(ggplot2::ggsave(paste0(base, ".png"), p, width = plot_base_width*2.0, height = plot_base_height*1.25, dpi = 400, limitsize=FALSE))
  suppressWarnings(ggplot2::ggsave(paste0(base, ".pdf"), p, width = plot_base_width*2.0, height = plot_base_height*1.25, limitsize=FALSE))
}

# --- Sankeys ---
plot_sankey_pop_series_withbars <- function(freq_mat_df, out_png, pop_order=NULL, label_threshold=0.01, bar_halfwidth=0.33, plot_title=NULL){
  long <- freq_mat_df %>% rename(group = 1) %>% pivot_longer(-group, names_to = "component", values_to = "value")
  if (is.null(pop_order)) pop_order <- order_pop_by_dominant_hap(freq_mat_df)
  long$group <- factor(long$group, levels = pop_order)
  long <- long %>% group_by(group) %>% mutate(value = value/sum(value)) %>% ungroup()
  components <- sort(unique(long$component))
  pal <- get_palette_for_levels(components, min_n = max(256, length(components)))
  stack_df <- long %>% arrange(group, component) %>% group_by(group) %>% mutate(ymin=cumsum(lag(value,default=0)), ymax=ymin+value) %>% ungroup()
  group_pos <- data.frame(group=levels(long$group), x=seq_along(levels(long$group)))
  stack_df <- left_join(stack_df, group_pos, by="group")
  rect_df <- stack_df %>% mutate(xmin=x-bar_halfwidth, xmax=x+bar_halfwidth, fill_col=pal[component])
  flows <- list()
  for(h in components){
    s <- rect_df %>% filter(component==h) %>% arrange(x); if(nrow(s)<2) next
    for(i in 1:(nrow(s)-1)){
      left  <- s[i, ]; right <- s[i+1, ]
      flows[[length(flows)+1]] <- data.frame(
        component = h,
        x = c(left$xmax, left$xmax, right$xmin, right$xmin),
        y = c(left$ymin, left$ymax, right$ymax, right$ymin),
        seg = paste0(h,"_",i),
        fill_col = pal[h]
      )
    }
  }
  flow_df <- if (length(flows)) do.call(rbind, flows) else data.frame(component=character(), x=numeric(), y=numeric(), seg=character(), fill_col=character())
  label_df <- rect_df %>% filter(value>=label_threshold) %>% mutate(ymid=(ymin+ymax)/2,label_txt=component)
  p <- ggplot()+
    geom_polygon(data=flow_df, aes(x=x,y=y,group=seg), fill=flow_df$fill_col, alpha=0.45, color=NA)+
    geom_rect(data=rect_df, aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax), fill=rect_df$fill_col, color="grey25", linewidth=0.25, alpha=0.98)+
    geom_text(data=label_df, aes(x=xmax - 0.01, y=ymid, label=label_txt), hjust=1, size=2.8, fontface="bold")+
    scale_x_continuous(breaks=group_pos$x,labels=group_pos$group,position="bottom")+
    scale_y_continuous(limits=c(0,1),expand=expansion(mult=0.02), labels=scales::percent_format(accuracy=1))+
    labs(x=NULL,y="Proportion", title=plot_title %||% "Population-series haplogroup flows (with bars)")+
    theme_minimal(base_size=13)+
    theme(panel.grid=element_blank(), axis.text.x=element_text(angle=60,hjust=1,vjust=1),
          panel.background=element_rect(fill="white",color=NA), plot.background=element_rect(fill="white",color=NA),
          plot.title=element_text(face="bold",size=12, hjust=0))+
    coord_cartesian(clip="off")
  base <- tools::file_path_sans_ext(out_png)
  ggsave(paste0(base,".png"),p,width=nature_sankey_width,height=nature_sankey_height,dpi=400, limitsize=FALSE)
  ggsave(paste0(base,".pdf"),p,width=nature_sankey_width,height=nature_sankey_height, limitsize=FALSE)
}
plot_sankey_pop_series_pure <- function(freq_mat_df, out_png, pop_order=NULL, label_threshold=0.01, bar_halfwidth=0.32, plot_title=NULL){
  long <- freq_mat_df %>% rename(group = 1) %>% pivot_longer(-group, names_to = "component", values_to = "value")
  if (is.null(pop_order)) pop_order <- order_pop_by_dominant_hap(freq_mat_df)
  long$group <- factor(long$group, levels = pop_order)
  long <- long %>% group_by(group) %>% mutate(value = value / sum(value)) %>% ungroup()
  components <- sort(unique(long$component))
  pal <- get_palette_for_levels(components, min_n = max(256, length(components)))
  stack_df <- long %>% arrange(group, component) %>% group_by(group) %>% mutate(ymin = cumsum(lag(value, default = 0)), ymax = ymin + value) %>% ungroup()
  group_pos <- data.frame(group = levels(long$group), x = seq_along(levels(long$group)), stringsAsFactors = FALSE)
  stack_df <- stack_df %>% left_join(group_pos, by = "group")
  flow_polys <- list(); label_rows <- list()
  for (comp in components) {
    sub_comp <- stack_df %>% filter(component == comp) %>% arrange(x); if (nrow(sub_comp) < 2) next
    for (i in 1:(nrow(sub_comp)-1)) {
      x0 <- sub_comp$x[i]; x1 <- sub_comp$x[i+1]
      ymin0 <- sub_comp$ymin[i]; ymax0 <- sub_comp$ymax[i]
      ymin1 <- sub_comp$ymin[i+1]; ymax1 <- sub_comp$ymax[i+1]
      v0 <- sub_comp$value[i];    v1 <- sub_comp$value[i+1]
      seg_id <- paste0(comp, "_", x0, "_", x1)
      poly_df <- data.frame(component = as.character(comp),
                            x = c(x0, x0, x1, x1),
                            y = c(ymin0, ymax0, ymax1, ymin1),
                            seg_id = seg_id, stringsAsFactors = FALSE)
      flow_polys[[length(flow_polys) + 1]] <- poly_df
      mid_x  <- (x0 + x1)/2; mid_y0 <- (ymin0 + ymax0)/2; mid_y1 <- (ymin1 + ymax1)/2; mid_y  <- (mid_y0 + mid_y1)/2
      label_rows[[length(label_rows) + 1]] <- data.frame(component = as.character(comp), seg_id = seg_id, x = mid_x, y = mid_y, v0 = v0, v1 = v1, stringsAsFactors = FALSE)
    }
  }
  flow_polys_df <- if (length(flow_polys)) do.call(rbind, flow_polys) else data.frame(component=character(), x=numeric(), y=numeric(), seg_id=character())
  flow_polys_df$fill_col <- pal[ flow_polys_df$component ]
  label_df <- if (length(label_rows)) do.call(rbind, label_rows) else data.frame(component=character(), seg_id=character(), x=numeric(), y=numeric(), v0=numeric(), v1=numeric())
  if (nrow(label_df)){
    label_df <- label_df %>% mutate(max_v = pmax(v0, v1), show  = max_v >= label_threshold, label_txt = component) %>% filter(show)
  }
  p <- ggplot() +
    geom_polygon(data = flow_polys_df, aes(x = x, y = y, group = seg_id), fill  = flow_polys_df$fill_col, alpha = 0.55, color = NA) +
    geom_text(data = label_df, aes(x = x, y = y, label = label_txt), size = 2.8, color = "black", fontface = "bold", lineheight = 0.9) +
    scale_x_continuous(breaks = group_pos$x, labels = group_pos$group, position = "bottom", expand = expansion(mult = 0.08)) +
    scale_y_continuous(limits = c(0,1), expand = expansion(mult = 0.02), labels = scales::percent_format(accuracy = 1)) +
    labs(x = NULL, y = "Proportion", title = plot_title %||% "Population-series haplogroup flows (pure)") +
    theme_minimal(base_size = 13) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
          plot.margin = margin(10, 20, 10, 20), panel.background = element_rect(fill = "white", color = NA),
          plot.background  = element_rect(fill = "white", color = NA), plot.title = element_text(face="bold", hjust=0))
  base <- tools::file_path_sans_ext(out_png)
  ggsave(paste0(base, ".png"), p, width = nature_sankey_width, height = nature_sankey_height, dpi = 400, limitsize=FALSE)
  ggsave(paste0(base, ".pdf"), p, width = nature_sankey_width, height = nature_sankey_height, limitsize=FALSE)
}

# --- Hap spread ---
compute_hap_spread_stats <- function(count_df, freq_df){
  base_ct <- count_df %>% group_by(hap) %>% summarise(total_n = sum(n), n_pop = n(), .groups = "drop")
  spread <- count_df %>% group_by(hap) %>% mutate(p = n/sum(n)) %>%
    summarise(shannon_pop = shannon(p), evenness_pop = ifelse(n() > 1, shannon_pop/log(n()), 1), .groups = "drop")
  pop_detail <- count_df %>% group_by(hap) %>% arrange(desc(n), population, .by_group=TRUE) %>%
    summarise(pop_detail = paste0(population, ":", n, collapse=";"),
              top_population = population[1], top_count = n[1], .groups = "drop")
  base_ct %>% left_join(spread, by="hap") %>% left_join(pop_detail, by="hap") %>% arrange(desc(total_n))
}

# --- Bootstrap helpers ---
get_pop_order_by_ORpeak <- function(df_boot) {
  df_boot <- df_boot %>% mutate(Haplogroup_up = toupper(Haplogroup), is_OR = grepl("^(O|R)", Haplogroup_up))
  pop_tot <- df_boot %>% count(population, name = "pop_total")
  or_cnt <- df_boot %>% filter(is_OR) %>% count(population, Haplogroup, name = "or_count") %>%
    left_join(pop_tot, by = "population") %>% mutate(or_freq = or_count / pop_total)
  best_or <- or_cnt %>% group_by(population) %>% summarise(max_or_freq = max(or_freq, na.rm = TRUE), best_or_hap = Haplogroup[which.max(or_freq)], .groups = "drop")
  pop_info <- pop_tot %>% left_join(best_or, by = "population") %>% mutate(max_or_freq = ifelse(is.na(max_or_freq), -1, max_or_freq)) %>%
    arrange(desc(max_or_freq), desc(pop_total), population)
  pop_info$population
}
bootstrap_once_with_div <- function(rep_id, df, n_total, n_per_sample) {
  idx <- sample.int(n_total, size = n_per_sample, replace = TRUE)
  sample_df <- df[idx, c("population", "Haplogroup"), drop = FALSE]
  counts <- sample_df %>% count(population, Haplogroup, name = "count")
  pop_counts <- sample_df %>% count(population, name = "pop_total")
  freq_tbl <- counts %>% left_join(pop_counts, by = "population") %>% mutate(freq = count / pop_total, rep  = rep_id) %>% select(rep, population, Haplogroup, freq)
  div_tbl <- freq_tbl %>% group_by(population, rep) %>% summarise(shannon = { p <- freq; p <- p[p > 0]; -sum(p * log(p)) }, simpson = { p <- freq; 1 - sum(p^2) }, .groups = "drop")
  list(freq_tbl = freq_tbl, div_tbl = div_tbl)
}
get_sig_table_for_hap <- function(subdat, hap_name) {
  pop_counts <- subdat %>% group_by(population) %>% summarise(n_points = n(), .groups="drop")
  pops_ok <- pop_counts %>% filter(n_points >= 2) %>% pull(population) %>% as.character()
  if (length(unique(pops_ok)) < 2) {
    return(tibble(Haplogroup=character(0), group1=character(0), group2=character(0), p=numeric(0), p.adj=numeric(0), p.adj.signif=character(0)))
  }
  sub_for_test <- subdat %>% filter(population %in% pops_ok)
  stat_test <- sub_for_test %>% rstatix::pairwise_wilcox_test(freq ~ population, p.adjust.method = "BH")
  if (nrow(stat_test) == 0) {
    return(tibble(Haplogroup=character(0), group1=character(0), group2=character(0), p=numeric(0), p.adj=numeric(0), p.adj.signif=character(0)))
  }
  stat_test %>% mutate(p.adj.signif = ifelse(is.na(p.adj.signif), "ns", p.adj.signif), Haplogroup = hap_name) %>%
    select(Haplogroup, group1, group2, p, p.adj, p.adj.signif)
}

# --- Plot single hap with optional significance ---
plot_hap_single <- function(subdat, pop_levels, pop_colors, n_per_sample, n_reps) {
  subdat <- subdat %>% group_by(population) %>% filter(any(freq > 0)) %>% ungroup()
  if (nrow(subdat) == 0) return(ggplot() + ggtitle("no data for this hap in any population"))
  used_pops <- intersect(pop_levels, unique(subdat$population))
  subdat$population <- factor(subdat$population, levels = used_pops, ordered = TRUE)
  y_lower <- min(subdat$freq, na.rm = TRUE); y_upper <- max(subdat$freq, na.rm = TRUE)
  base_plot <- ggplot(subdat, aes(x = population, y = freq)) +
    ggbeeswarm::geom_quasirandom(aes(color = population),
                                 method = qr_method, width = qr_width, alpha = qr_alpha, size = qr_size, show.legend = FALSE) +
    geom_boxplot(aes(color = population), width = 0.5, outlier.shape = NA, fill = NA, size = 0.6, show.legend = FALSE) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 0, color = "red", size = 0.6) +
    scale_color_manual(values = pop_colors[used_pops]) +
    labs(x = "Population", y = "Frequency", title = paste0(unique(subdat$Haplogroup), " (n=", n_per_sample, " /rep, ", n_reps, " reps)")) +
    theme_bw() + theme(axis.text.x  = element_text(angle = 30, hjust = 1, vjust = 1, size = 8),
                       axis.text.y  = element_text(size = 8),
                       axis.title.x = element_text(size = 9),
                       axis.title.y = element_text(size = 9),
                       plot.title   = element_text(hjust = 0.5, size = 9)) +
    coord_cartesian(ylim = c(y_lower, y_upper))

  if (!isTRUE(freqplot_show_sig)) return(base_plot)

  # Compute pairwise Wilcoxon (BH) on-the-fly
  stat_tbl <- tryCatch(
    subdat %>% rstatix::pairwise_wilcox_test(freq ~ population, p.adjust.method = "BH"),
    error = function(e) NULL
  )
  if (is.null(stat_tbl) || nrow(stat_tbl) == 0) return(base_plot)

  # Filter to desired comparisons
  if (!is.null(freqplot_sig_ref_pop) && freqplot_sig_ref_pop %in% levels(subdat$population)) {
    stat_tbl <- stat_tbl %>% filter(group1 == freqplot_sig_ref_pop | group2 == freqplot_sig_ref_pop)
  }
  stat_tbl <- stat_tbl %>% filter(!is.na(p.adj) & p.adj <= freqplot_sig_alpha) %>% arrange(p.adj)

  if (nrow(stat_tbl) == 0) return(base_plot)

  # Limit number of brackets
  if (!is.null(freqplot_sig_max_pairs) && nrow(stat_tbl) > freqplot_sig_max_pairs) {
    stat_tbl <- stat_tbl %>% slice_head(n = freqplot_sig_max_pairs)
  }

  comps <- lapply(seq_len(nrow(stat_tbl)), function(i) c(as.character(stat_tbl$group1[i]), as.character(stat_tbl$group2[i])))
  ymax <- max(subdat$freq, na.rm = TRUE)
  step <- max(0.06 * ymax, 0.02)
  y_positions <- seq(ymax + step, by = step, length.out = length(comps))
  labels <- ifelse(is.na(stat_tbl$p.adj.signif), "ns", stat_tbl$p.adj.signif)

  base_plot + ggsignif::geom_signif(comparisons = comps,
                                    annotations = labels,
                                    y_position = y_positions,
                                    tip_length = 0.01, vjust = 0.5, textsize = 3)
}

# --- Diversity scatter (quasirandom) ---
plot_diversity_metric <- function(diversity_long, metric = c("shannon","simpson"), pop_levels, pop_colors) {
  metric <- match.arg(metric)
  pop_order_by_mean <- diversity_long %>% group_by(population) %>% summarise(mv = mean(.data[[metric]], na.rm = TRUE), .groups = "drop") %>% arrange(desc(mv)) %>% pull(population)
  diversity_long$population <- factor(diversity_long$population, levels = pop_order_by_mean, ordered = TRUE)
  y_lower <- min(diversity_long[[metric]], na.rm = TRUE); y_upper <- max(diversity_long[[metric]], na.rm = TRUE)
  p <- ggplot(diversity_long, aes(x = population, y = .data[[metric]])) +
    ggbeeswarm::geom_quasirandom(aes(color = population),
                                 method = qr_method, width = qr_width, alpha = qr_alpha, size = qr_size, show.legend = FALSE) +
    geom_boxplot(aes(color = population), width = 0.5, outlier.shape = NA, fill = NA, size = 0.6, show.legend = FALSE) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 0, color = "red", size = 0.6) +
    scale_color_manual(values = pop_colors[pop_order_by_mean]) +
    labs(x = "Population", y = metric, title = paste0("Population diversity (", metric, ")")) +
    theme_bw() + theme(axis.text.x  = element_text(angle = 30, hjust = 1, vjust = 1, size = 8),
                       axis.text.y  = element_text(size = 8),
                       axis.title.x = element_text(size = 9),
                       axis.title.y = element_text(size = 9),
                       plot.title   = element_text(hjust = 0.5, size = 10)) +
    coord_cartesian(ylim = c(y_lower, y_upper))
  p
}

# --- Facets ---
plot_and_save_facets <- function(df_in, suffix_label, make_title, outdir_facets, pop_levels, nature_palette, ncol_facets, per_col_w, per_row_h) {
  if (nrow(df_in) == 0) { message("facet skipped: ", suffix_label, " has no data"); return(invisible(NULL)) }
  df_in <- df_in %>% mutate(population  = factor(population, levels = pop_levels, ordered = TRUE),
                            HaplogroupF = factor(Haplogroup, ordered = FALSE))
  median_by_pop <- df_in %>% group_by(HaplogroupF, population) %>% summarise(freq_median = median(freq, na.rm = TRUE), .groups = "drop") %>% mutate(pop_index = as.numeric(factor(population, levels = pop_levels, ordered = TRUE)))
  corr_table <- median_by_pop %>% group_by(HaplogroupF) %>%
    summarise(cor_trend = ifelse(n() >= 2 && sd(pop_index) > 0 && sd(freq_median) > 0, cor(pop_index, freq_median, method = "pearson"), NA_real_), .groups = "drop") %>%
    arrange(is.na(cor_trend), cor_trend)
  hap_levels_sorted <- as.character(corr_table$HaplogroupF)
  df_in <- df_in %>% mutate(HaplogroupF = factor(HaplogroupF, levels = hap_levels_sorted, ordered = TRUE))
  pop_colors <- setNames(rep(nature_palette, length.out = length(pop_levels)), pop_levels)

  p_facets <- ggplot(df_in, aes(x = population, y = freq)) +
    ggbeeswarm::geom_quasirandom(aes(color = population),
                                 method = qr_method, width = qr_width, alpha = qr_alpha, size = qr_size, show.legend = FALSE) +
    geom_boxplot(aes(color = population), width = 0.5, outlier.shape = NA, fill = NA, size = 0.4, show.legend = FALSE) +
    stat_summary(fun = median, geom = "crossbar", width = 0.5, fatten = 0, color = "red", size = 0.4) +
    scale_color_manual(values = pop_colors) +
    labs(x = "Population", y = "Frequency", title = make_title) +
    theme_bw() +
    theme(strip.text = element_text(size = 6, face = "bold"),
          axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1, size = 5),
          axis.text.y = element_text(size = 5),
          axis.title.x = element_text(size = 7),
          axis.title.y = element_text(size = 7),
          plot.title   = element_text(hjust = 0.5, size = 10),
          panel.spacing.x = unit(0.4, "lines"),
          panel.spacing.y = unit(0.4, "lines")) +
    facet_wrap(~ HaplogroupF, ncol = ncol_facets, scales = "free_y")

  print(p_facets)
  n_haps <- length(unique(df_in$HaplogroupF)); n_rows <- ceiling(n_haps / ncol_facets)
  fig_w  <- per_col_w * ncol_facets; fig_h  <- per_row_h * n_rows
  png_path <- file.path(outdir_facets, paste0("facet_", suffix_label, ".png"))
  pdf_path <- file.path(outdir_facets, paste0("facet_", suffix_label, ".pdf"))
  ggsave(png_path, plot = p_facets, width = fig_w, height = fig_h, dpi = 300, limitsize=FALSE)
  ggsave(pdf_path, plot = p_facets, width = fig_w, height = fig_h, limitsize=FALSE)

  corr_outfile <- file.path(outdir_facets, paste0("haplogroup_correlation_order_", suffix_label, ".csv"))
  readr::write_csv(corr_table %>% rename(Haplogroup = HaplogroupF), corr_outfile)
}

# --- Geo tests ---
run_optional_geo_tests <- function(D, meta, out_prefix){
  if (!is.data.frame(meta) || is.null(colnames(meta)) || ncol(meta)==0) return(invisible(NULL))
  D <- sanitize_dist_matrix(D, name = out_prefix)
  colnames(meta) <- tolower(colnames(meta))
  meta <- meta %>% distinct(population, .keep_all=TRUE) %>% filter(population %in% rownames(D))
  if (nrow(meta) < 3) return(invisible(NULL))
  if (all(c("lat","lon") %in% colnames(meta))){
    coords <- meta[match(rownames(D), meta$population), c("lat","lon")]
    G <- as.matrix(dist(coords))
    mantel_res <- vegan::mantel(as.dist(D), as.dist(G), permutations=9999)
    writeLines(c(paste0("Mantel r = ", signif(mantel_res$statistic,4)),
                 paste0("Mantel p = ", signif(mantel_res$signif,4))), con=paste0(out_prefix, "_Mantel.txt"))
  }
  if ("region" %in% colnames(meta)){
    df <- data.frame(population = rownames(D)) %>% left_join(meta, by="population")
    ad <- vegan::adonis2(as.dist(D) ~ region, data=df, permutations=9999)
    capture.output(ad, file=paste0(out_prefix, "_PERMANOVA.txt"))
  }
}

# --- Bootstrap runner ---
run_bootstrap_for_hh <- function(df_hh, HH, out_dir,
                                 n_per_sample   = bootstrap_n_per_sample,
                                 n_reps         = bootstrap_n_reps,
                                 pop_levels_opt = bootstrap_pop_levels,
                                 ncol_facets    = bootstrap_ncol_facets,
                                 per_col_w      = bootstrap_per_col_w,
                                 per_row_h      = bootstrap_per_row_h,
                                 keep_only      = bootstrap_keep_only_levels) {
  message2("Start Bootstrap (HH=%d)", HH)
  outdir_main   <- file.path(out_dir, "bootstrap")
  outdir_facets <- file.path(outdir_main, "facet_pages")
  dir.create(outdir_main,   showWarnings = FALSE, recursive = TRUE)
  dir.create(outdir_facets, showWarnings = FALSE, recursive = TRUE)

  df_boot <- df_hh %>% transmute(id = sample_id %||% NA, Haplogroup_full = haplogroup_full %||% NA,
                                 population = population, Haplogroup = hap) %>% filter(!is.na(Haplogroup), !is.na(population))

  if (is.null(pop_levels_opt) || length(pop_levels_opt) == 0) pop_levels <- get_pop_order_by_ORpeak(df_boot) else pop_levels <- pop_levels_opt
  if (keep_only) df_boot <- df_boot %>% filter(population %in% pop_levels)

  n_total <- nrow(df_boot); if (n_total == 0) { message2("HH=%d has no rows for bootstrap, skip.", HH); return(invisible(NULL)) }
  message2("HH=%d bootstrap total rows: %d", HH, n_total)

  nature_palette <- c("#4C72B0", "#55A868", "#C44E52", "#8172B2", "#CCB974", "#64B5CD", "#8C8C8C", "#E17C05", "#937860", "#5A9E6F")
  pop_colors <- setNames(rep(nature_palette, length.out = length(pop_levels)), pop_levels)

  boot_list <- vector("list", n_reps)
  for (r in seq_len(n_reps)) boot_list[[r]] <- bootstrap_once_with_div(r, df_boot, n_total, n_per_sample)

  freq_long      <- bind_rows(lapply(boot_list, `[[`, "freq_tbl"))
  diversity_long <- bind_rows(lapply(boot_list, `[[`, "div_tbl"))

  freq_table_csv      <- file.path(outdir_main, "haplogroup_bootstrap_freq_long.csv")
  sig_table_csv       <- file.path(outdir_main, "haplogroup_significance_results.csv")
  diversity_table_csv <- file.path(outdir_main, "population_diversity_bootstrap_long.csv")
  write.csv(freq_long,      freq_table_csv,      row.names = FALSE)
  write.csv(diversity_long, diversity_table_csv, row.names = FALSE)

  hap_levels_thisHH <- sort(unique(freq_long$Haplogroup))
  all_sig_rows <- list()
  for (hap in hap_levels_thisHH) {
    subdat <- freq_long %>% filter(Haplogroup == hap); if (nrow(subdat) == 0) next
    sig_tab <- get_sig_table_for_hap(subdat, hap); all_sig_rows[[hap]] <- sig_tab

    p_hap <- plot_hap_single(subdat=subdat, pop_levels=pop_levels, pop_colors=pop_colors,
                             n_per_sample=n_per_sample, n_reps=n_reps)
    w <- max(3.5, 0.4 * length(pop_levels)); h <- 3.5
    ggsave(file.path(outdir_main, paste0("freqplot_", hap, ".png")), plot = p_hap, width = w, height = h, dpi = 300, limitsize=FALSE)
    ggsave(file.path(outdir_main, paste0("freqplot_", hap, ".pdf")), plot = p_hap, width = w, height = h, limitsize=FALSE)
  }
  sig_results_all <- bind_rows(all_sig_rows); write.csv(sig_results_all, sig_table_csv, row.names = FALSE)

  p_shannon <- plot_diversity_metric(diversity_long, metric="shannon", pop_levels=pop_levels, pop_colors=pop_colors)
  p_simpson <- plot_diversity_metric(diversity_long, metric="simpson",  pop_levels=pop_levels, pop_colors=pop_colors)
  ggsave(file.path(outdir_main, "diversity_shannon.png"), plot = p_shannon, width = 6, height = 4, dpi = 300, limitsize=FALSE)
  ggsave(file.path(outdir_main, "diversity_shannon.pdf"),  plot = p_shannon, width = 6, height = 4, limitsize=FALSE)
  ggsave(file.path(outdir_main, "diversity_simpson.png"),  plot = p_simpson, width = 6, height = 4, dpi = 300, limitsize=FALSE)
  ggsave(file.path(outdir_main, "diversity_simpson.pdf"),  plot = p_simpson, width = 6, height = 4, limitsize=FALSE)

  freq_long_forFacet <- freq_long %>% mutate(Haplogroup = as.character(Haplogroup), hap_len = nchar(Haplogroup),
                                             population = factor(population, levels = pop_levels, ordered = TRUE))
  freq_lenK <- freq_long_forFacet %>% filter(hap_len == HH)
  plot_and_save_facets(df_in=freq_lenK, suffix_label=paste0("length", sprintf("%02d", HH), "_sorted"),
                       make_title=paste0("Haplogroup bootstrap frequency (length=", HH, ", sorted)"),
                       outdir_facets=outdir_facets, pop_levels=pop_levels, nature_palette=nature_palette,
                       ncol_facets=bootstrap_ncol_facets, per_col_w=bootstrap_per_col_w, per_row_h=bootstrap_per_row_h)

  n_pops_total <- length(pop_levels)
  haps_with_all <- freq_lenK %>% group_by(Haplogroup) %>% summarise(n_pop = n_distinct(population), .groups = "drop") %>%
    filter(n_pop == n_pops_total) %>% pull(Haplogroup)
  freq_lenK_all <- freq_lenK %>% filter(Haplogroup %in% haps_with_all)
  plot_and_save_facets(df_in=freq_lenK_all, suffix_label=paste0("length", sprintf("%02d", HH), "_sorted_allPops"),
                       make_title=paste0("Haplogroup bootstrap frequency (length=", HH, ", present in all ", n_pops_total, " groups)"),
                       outdir_facets=outdir_facets, pop_levels=pop_levels, nature_palette=nature_palette,
                       ncol_facets=bootstrap_ncol_facets, per_col_w=bootstrap_per_col_w, per_row_h=bootstrap_per_row_h)
  message2("Bootstrap done (HH=%d), output: %s", HH, normalizePath(outdir_main))
}

# --- Main ---
message2("Reading input: %s", input_xlsx)
df0 <- safe_read_input(input_xlsx)
df0$population <- sanitize_label_ascii(df0$population)

# NEW: Try to load haplogroup hierarchy from file; if not found, fall back to character truncation
df_hap_hierarchy <- read_haplogroup_hierarchy(haplogroup_hierarchy_csv)

if (!is.null(df_hap_hierarchy)) {
  # Use hierarchy from file
  message2("Using haplogroup hierarchy from file: %s", haplogroup_hierarchy_csv)
  df0 <- merge_haplogroup_hierarchy(df0, df_hap_hierarchy)
  # Determine HH_range from the number of hap_level_* columns available
  hap_level_cols <- grep("^hap_level_", colnames(df0), value = TRUE)
  if (length(hap_level_cols) > 0) {
    HH_range <- seq_len(length(hap_level_cols))
    message2("Auto-detected %d haplogroup levels from file; HH_range set to: %s", length(hap_level_cols), paste(HH_range, collapse=", "))
  }
} else {
  # Fall back: use character truncation from haplogroup_full (legacy mode)
  message2("Haplogroup hierarchy file not found or invalid; falling back to character truncation mode")
  if ("haplogroup_full" %in% colnames(df0)) {
    df0$hap_clean <- sanitize_hap(df0$haplogroup_full)
    message2("Generated hap_level_* columns via character truncation")
    # Create hap_level columns from truncation
    for (lv in HH_range) {
      df0[[paste0("hap_level_", lv)]] <- hap_substr(df0$hap_clean, lv)
    }
  } else {
    stop("Neither haplogroup_hierarchy_csv file nor haplogroup_full column found in input!")
  }
}

if (use_equalize){
  pop_sizes <- table(df0$population); target_n <- min(pop_sizes)
  message2("Equalize sampling enabled; target n per population = %d", target_n); set.seed(1)
  df_in <- do.call(rbind, lapply(split(df0, df0$population), function(tab){
    n <- nrow(tab); if (n > target_n) tab[sample.int(n, target_n), , drop=FALSE] else tab
  }))
} else df_in <- df0

# Verify that hap_level columns exist
hap_level_cols <- grep("^hap_level_", colnames(df_in), value = TRUE)
if (length(hap_level_cols) == 0) stop("No hap_level_* columns found after loading; check input.")
HH_range <- seq_len(length(hap_level_cols))  # Ensure HH_range matches available columns

# meta (optional)
meta_in <- NULL
if (file.exists(meta_csv)){
  meta_in <- tryCatch(read.csv(meta_csv, stringsAsFactors=FALSE, fileEncoding="UTF-8"), error=function(e) NULL)
  if (!is.null(meta_in)){
    meta_in$population <- sanitize_label_ascii(meta_in$population)
    if ("group" %in% names(meta_in))  meta_in$group  <- sanitize_label_ascii(meta_in$group)
    if ("region" %in% names(meta_in)) meta_in$region <- sanitize_label_ascii(meta_in$region)
    message2("Loaded meta: %s", meta_csv)
  }
}

# Y-STR (optional)
if (file.exists(ystr_csv)){
  message2("Y-STR table detected: %s, computing Rst", ystr_csv)
  ystr_df <- tryCatch(read.csv(ystr_csv, stringsAsFactors=FALSE, fileEncoding="UTF-8"), error=function(e) NULL)
  if (!is.null(ystr_df)){
    ystr_df$population <- sanitize_label_ascii(ystr_df$population)
    out_dir_rst <- file.path(root_out, "out_Rst"); dir.create(out_dir_rst, showWarnings=FALSE, recursive=TRUE)
    D_RST <- compute_pairwise_Rst(ystr_df, pop_col="population")
    write.csv(D_RST, file.path(out_dir_rst, "pairwise_Rst.csv"), fileEncoding="UTF-8")
    plot_distance_heatmap(D_RST, file.path(out_dir_rst, "heatmap_Rst.png"), title="Rst heatmap (Y-STR)")
    plot_dendrogram(D_RST, file.path(out_dir_rst, "dendrogram_Rst.pdf"))
    plot_nj_tree(D_RST, file.path(out_dir_rst, "NJ_Rst.pdf"), title="NJ tree (Rst, Y-STR)")
    run_mds_plot(D_RST, file.path(out_dir_rst, "MDS_Rst.png"))
  }
}

for (HH in HH_range){
  message2("==== Processing HH=%s ====", HH)
  out_dir <- file.path(root_out, sprintf("out_HH%02d", HH))
  dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
  dir.create(file.path(out_dir,"tables"), showWarnings=FALSE, recursive=TRUE)
  dir.create(file.path(out_dir,"figs"),   showWarnings=FALSE, recursive=TRUE)

  # NEW: Use hap_level_* column directly instead of truncating
  hap_col_name <- paste0("hap_level_", HH)
  if (!hap_col_name %in% colnames(df_in)) {
    message2("WARNING: Column %s not found; skipping HH=%d", hap_col_name, HH)
    next
  }
  df <- df_in %>% mutate(hap = .data[[hap_col_name]]) %>% filter(!is.na(hap) & nchar(hap) > 0)

  # freq/count matrices
  freq_df     <- freq_table(df, hap_col="hap", by_col="population")
  freq_mat_df <- make_freq_matrix(freq_df, by_col="population", hap_col="hap", val_col="freq")
  count_df    <- df %>% count(population, hap, name="n")
  count_mat_df<- make_freq_matrix(count_df, by_col="population", hap_col="hap", val_col="n")

  write.csv(freq_mat_df,  file.path(out_dir,"tables",sprintf("freq_matrix_HH%02d.csv",  HH)), row.names=FALSE, fileEncoding="UTF-8")
  write.csv(count_mat_df, file.path(out_dir,"tables",sprintf("count_matrix_HH%02d.csv", HH)), row.names=FALSE, fileEncoding="UTF-8")

  # hap spread
  hap_spread_tab <- compute_hap_spread_stats(count_df, freq_df)
  write.csv(hap_spread_tab, file.path(out_dir,"tables",sprintf("hap_subclade_diversity_HH%02d.csv", HH)), row.names=FALSE, fileEncoding="UTF-8")

  # bipartite edges
  bip_edges <- count_df %>% rename(count = n) %>% left_join(freq_df %>% select(population, hap, freq), by=c("population","hap")) %>% arrange(desc(freq))
  write.csv(bip_edges, file.path(out_dir,"tables",sprintf("bipartite_edges_HH%02d.csv", HH)), row.names=FALSE, fileEncoding="UTF-8")

  # diversity indices
  div_list <- lapply(split(freq_df, freq_df$population), function(x){
    p <- x$freq
    data.frame(population = x$population[1], k = nrow(x),
               shannon = shannon(p), simpson = simpson(p),
               heff_expH = heff_shannon(p), heff_1sumP2 = heff_simpson(p))
  })
  div_tab <- do.call(rbind, div_list)
  write.csv(div_tab, file.path(out_dir,"tables",sprintf("diversity_indices_HH%02d.csv", HH)), row.names=FALSE, fileEncoding="UTF-8")

  # prevalence
  classify_prevalence <- function(freq_df, count_df){
    x <- left_join(freq_df, count_df, by=c("population","hap"))
    x$class <- with(x, ifelse(n.y == 1, "Singleton",
                              ifelse(n.y == 2, "Doubleton",
                                     ifelse(freq <= thr_rare, "Very rare",
                                            ifelse(freq <= thr_common, "Rare", "Common")))))
    detail <- x %>% select(population, hap, count = n.y, freq, class)
    summary <- detail %>% count(population, class, name="hap_n") %>% pivot_wider(names_from=class, values_from=hap_n, values_fill=0)
    list(detail=detail, summary=summary)
  }
  pv <- classify_prevalence(freq_df, count_df)
  write.csv(pv$detail,  file.path(out_dir,"tables",sprintf("hap_prevalence_classes_HH%02d.csv", HH)), row.names=FALSE, fileEncoding="UTF-8")
  write.csv(pv$summary, file.path(out_dir,"tables",sprintf("hap_prevalence_summary_by_population_HH%02d.csv", HH)), row.names=FALSE, fileEncoding="UTF-8")

  # enrichment
  enr <- hap_enrichment(df, hap_col="hap", pop_col="population")
  write.csv(enr, file.path(out_dir,"tables",sprintf("enrichment_Fisher_HH%02d.csv", HH)), row.names=FALSE, fileEncoding="UTF-8")

  # distances
  mat   <- as.matrix(freq_mat_df[,-1]); rownames(mat) <- freq_mat_df$population
  D_TVD <- pairwise_matrix(mat, TVD)
  D_JSD <- pairwise_matrix(mat, function(a,b) JSD_pair(a,b))
  D_NEI <- pairwise_neiD(mat)
  D_DA  <- pairwise_DA(mat)
  write.csv(D_TVD, file.path(out_dir,"tables",sprintf("pairwise_TVD_HH%02d.csv", HH)), fileEncoding="UTF-8")
  write.csv(D_JSD, file.path(out_dir,"tables",sprintf("pairwise_JSDist_HH%02d.csv", HH)), fileEncoding="UTF-8")
  write.csv(D_NEI, file.path(out_dir,"tables",sprintf("pairwise_NeiD_HH%02d.csv", HH)), fileEncoding="UTF-8")
  write.csv(D_DA,  file.path(out_dir,"tables",sprintf("pairwise_DA_HH%02d.csv",  HH)), fileEncoding="UTF-8")

  # heatmap + trees + NJ
  plot_distance_heatmap(D_TVD, file.path(out_dir,"figs",sprintf("heatmap_TVD_HH%02d.png", HH)), title=sprintf("TVD heatmap (HH=%d)", HH))
  plot_distance_heatmap(D_JSD, file.path(out_dir,"figs",sprintf("heatmap_JSD_HH%02d.png", HH)), title=sprintf("JSDist heatmap (HH=%d)", HH))
  plot_distance_heatmap(D_NEI, file.path(out_dir,"figs",sprintf("heatmap_NeiD_HH%02d.png", HH)), title=sprintf("Nei's D heatmap (HH=%d)", HH))
  plot_distance_heatmap(D_DA,  file.path(out_dir,"figs",sprintf("heatmap_DA_HH%02d.png",  HH)), title=sprintf("DA heatmap (HH=%d)", HH))
  plot_dendrogram(D_TVD, file.path(out_dir,"figs",sprintf("dendrogram_TVD_HH%02d.pdf", HH)))
  plot_dendrogram(D_NEI, file.path(out_dir,"figs",sprintf("dendrogram_NeiD_HH%02d.pdf", HH)))
  plot_dendrogram(D_DA,  file.path(out_dir,"figs",sprintf("dendrogram_DA_HH%02d.pdf",  HH)))
  plot_nj_tree(D_NEI, file.path(out_dir,"figs",sprintf("NJ_NeiD_HH%02d.pdf", HH)), title=sprintf("NJ tree (Nei's D, HH=%d)", HH))
  plot_nj_tree(D_DA,  file.path(out_dir,"figs",sprintf("NJ_DA_HH%02d.pdf",  HH)),  title=sprintf("NJ tree (DA, HH=%d)", HH))

  # MDS
  run_mds_plot(D_TVD, file.path(out_dir,"figs",sprintf("MDS_TVD_HH%02d.png", HH)))
  run_mds_plot(D_JSD, file.path(out_dir,"figs",sprintf("MDS_JSD_HH%02d.png", HH)))
  run_mds_plot(D_NEI, file.path(out_dir,"figs",sprintf("MDS_NeiD_HH%02d.png", HH)))
  run_mds_plot(D_DA,  file.path(out_dir,"figs",sprintf("MDS_DA_HH%02d.png",  HH)))

  # PCA / CA
  run_pca_plot(freq_mat_df, file.path(out_dir,"figs",sprintf("PCA_HH%02d.png", HH)),
               file.path(out_dir,"tables",sprintf("PCA_scores_HH%02d.csv", HH)))
  run_ca_plot(count_mat_df, file.path(out_dir,"figs",sprintf("CA_HH%02d.png",  HH)),
              file.path(out_dir,"tables",sprintf("CA_coords_HH%02d.csv", HH)))

  # stacked bars
  plot_stacked_bars(freq_mat_df, file.path(out_dir,"figs",sprintf("StackedBars_HH%02d.png", HH)))

  # population order
  pop_order_opt <- order_pop_by_dominant_hap(freq_mat_df)

  # Two Sankeys
  sankey_pure_file <- file.path(out_dir,"figs",sprintf("Sankey_pure_HH%02d.png", HH))
  plot_sankey_pop_series_pure(freq_mat_df, out_png=sankey_pure_file, pop_order=pop_order_opt, label_threshold=0.01, bar_halfwidth=0.32,
                              plot_title=sprintf("Haplogroup flows across populations (pure, HH=%d)", HH))
  sankey_bar_file <- file.path(out_dir,"figs",sprintf("Sankey_withBars_HH%02d.png", HH))
  plot_sankey_pop_series_withbars(freq_mat_df, out_png=sankey_bar_file, pop_order=pop_order_opt, label_threshold=0.01, bar_halfwidth=0.33,
                                  plot_title=sprintf("Haplogroup flows across populations (with bars, HH=%d)", HH))

  # single-layer alluvial
  plot_alluvial_simple(df, hap_col="hap", pop_col="population", out_png=file.path(out_dir,"figs",sprintf("Alluvial_simple_HH%02d.png", HH)))
  plot_alluvial_simple_freq(freq_mat_df, out_png=file.path(out_dir,"figs",sprintf("Alluvial_simpleFreq_HH%02d.png", HH)))

  # multilayer alluvial
  if (make_multilayer_sankey && HH >= 2){
    HC <- sort(unique(pmin(seq_len(HH), multilayer_max_k)))
    mdat <- make_multilayer_long(df, HH_chain=HC)

    plot_alluvial_multilayer_counts_from_long(mdat, out_png=file.path(out_dir,"figs", sprintf("Alluvial_multilayer_chain_toHH%02d.png", HH)))
    plot_alluvial_multilayer_freq_from_long(mdat,   out_png=file.path(out_dir,"figs", sprintf("Alluvial_multilayerFreq_chain_toHH%02d.png", HH)))

    root_freq_tab <- mdat$long %>% group_by(.data[[mdat$axes[1]]]) %>% summarise(global_freq = sum(freq), .groups="drop") %>% arrange(desc(global_freq))
    root_to_plot <- root_freq_tab %>% filter(global_freq >= perhap_min_freq) %>% slice_head(n = perhap_max_plot) %>% pull(1)
    if (length(root_to_plot)){
      for (rh in root_to_plot){
        out_raw <- file.path(out_dir,"figs", sprintf("Alluvial_multilayerFreq_chain_toHH%02d_root_%s_raw.png", HH, rh))
        plot_alluvial_multilayer_freq_singlehap_raw(mdat, rh, out_raw)
        out_scaled <- file.path(out_dir,"figs", sprintf("Alluvial_multilayerFreq_chain_toHH%02d_root_%s_scaled.png", HH, rh))
        plot_alluvial_multilayer_freq_singlehap_scaled(mdat, rh, out_scaled)
      }
    }
  }

  # geo
  if (!is.null(meta_in)){
    try(run_optional_geo_tests(D_TVD, meta_in,
                               out_prefix=file.path(out_dir,"tables", sprintf("GeoTests_TVD_HH%02d", HH))), silent=TRUE)
  }

  # bootstrap
  if (bootstrap_enable) {
    run_bootstrap_for_hh(df_hh=df, HH=HH, out_dir=out_dir,
                         n_per_sample=bootstrap_n_per_sample, n_reps=bootstrap_n_reps,
                         pop_levels_opt=bootstrap_pop_levels, ncol_facets=bootstrap_ncol_facets,
                         per_col_w=bootstrap_per_col_w, per_row_h=bootstrap_per_row_h,
                         keep_only=bootstrap_keep_only_levels)
  }

  message2("HH=%d done, outputs at: %s", HH, normalizePath(out_dir))
}

if (!exists("release_all_handles")) {
  release_all_handles <- function(move_wd = FALSE) { graphics.off(); if (move_wd) { try(setwd(".."), silent = TRUE) } }
}
release_all_handles(move_wd = FALSE)
message2("All HH complete. Output dir: %s", normalizePath(root_out))
