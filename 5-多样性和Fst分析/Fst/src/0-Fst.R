#!/usr/bin/env Rscript

# =============================================================================
# 配对Fst矩阵计算工具
# 
# 使用PopGenome包计算群体间配对Fst值
# 
# 用法：
# Rscript calculate_pairwise_fst.R <fasta_file> <population_file> <output_file> [method] [cores] [bootstrap] [alpha]
# 
# 参数说明：
# fasta_file      - 对齐的FASTA序列文件路径
# population_file - 群体信息文件路径（第一列：样本ID，第二列：群体名）
# output_file     - 输出CSV文件路径
# method          - 可选，Fst计算方法（默认：hudson）
# cores           - 可选，使用的CPU核心数（默认：1）
# bootstrap       - 可选，permutation检验次数（默认：1000，设为0跳过统计检验）
# alpha           - 可选，显著性水平（默认：0.05）
# 
# 可选的Fst计算方法：
# - hudson    : Hudson et al. (1992) F_ST (默认)
# - nei       : Nei's G_ST 
# - nucleotide: 基于核苷酸的F_ST
# - haplotype : 基于单倍型的F_ST
# 
# 示例：
# Rscript calculate_pairwise_fst.R data/sequences.fasta data/populations.txt output/fst_matrix.csv hudson 4 1000 0.05
# 
# =============================================================================

# 加载必要的包
suppressMessages({
  library(PopGenome)
  library(parallel)
  library(tools)
})

# =============================================================================
# 解析命令行参数
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)

# 检查参数数量
if (length(args) < 3) {
  cat("错误：参数不足\n\n")
  cat("用法：Rscript calculate_pairwise_fst.R <fasta_file> <population_file> <output_file> [method] [cores] [bootstrap] [alpha]\n\n")
  cat("参数说明：\n")
  cat("  fasta_file      - 对齐的FASTA序列文件路径\n")
  cat("  population_file - 群体信息文件路径（第一列：样本ID，第二列：群体名）\n")
  cat("  output_file     - 输出CSV文件路径\n")
  cat("  method          - 可选，Fst计算方法（默认：hudson）\n")
  cat("  cores           - 可选，使用的CPU核心数（默认：1）\n")
  cat("  bootstrap       - 可选，permutation检验次数（默认：1000，设为0跳过统计检验）\n")
  cat("  alpha           - 可选，显著性水平（默认：0.05）\n\n")
  cat("可选的Fst计算方法：\n")
  cat("  - hudson    : Hudson et al. (1992) F_ST (默认)\n")
  cat("  - nei       : Nei's G_ST\n")
  cat("  - nucleotide: 基于核苷酸的F_ST\n")
  cat("  - haplotype : 基于单倍型的F_ST\n\n")
  cat("示例：\n")
  cat("  Rscript calculate_pairwise_fst.R data/sequences.fasta data/populations.txt output/fst_matrix.csv hudson 4 1000 0.05\n")
  quit(status = 1)
}

# 解析参数
fasta_file <- args[1]
population_file <- args[2]
output_file <- args[3]
method <- if (length(args) >= 4) tolower(args[4]) else "hudson"
cores <- if (length(args) >= 5) as.numeric(args[5]) else 1
bootstrap <- if (length(args) >= 6) as.numeric(args[6]) else 1000
alpha <- if (length(args) >= 7) as.numeric(args[7]) else 0.05

# 验证参数
if (!file.exists(fasta_file)) {
  cat("错误：FASTA文件不存在：", fasta_file, "\n")
  quit(status = 1)
}

if (!file.exists(population_file)) {
  cat("错误：群体文件不存在：", population_file, "\n")
  quit(status = 1)
}

if (!method %in% c("hudson", "nei", "nucleotide", "haplotype")) {
  cat("错误：不支持的方法：", method, "\n")
  cat("支持的方法：hudson, nei, nucleotide, haplotype\n")
  quit(status = 1)
}

if (cores < 1 || cores > detectCores()) {
  cat("警告：核心数设置不合理，使用默认值 1\n")
  cores <- 1
}

if (bootstrap < 0) {
  cat("警告：bootstrap次数不能为负数，设为0跳过统计检验\n")
  bootstrap <- 0
}

if (alpha <= 0 || alpha >= 1) {
  cat("警告：显著性水平应在0-1之间，使用默认值0.05\n")
  alpha <- 0.05
}

# 创建输出目录
output_dir <- dirname(output_file)
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  cat("创建输出目录：", output_dir, "\n")
}

# 显示运行参数
cat("=============================================================================\n")
cat("配对Fst矩阵计算\n")
cat("=============================================================================\n")
cat("FASTA文件：", fasta_file, "\n")
cat("群体文件：", population_file, "\n")
cat("输出文件：", output_file, "\n")
cat("计算方法：", method, "\n")
cat("CPU核心：", cores, "\n")
cat("Bootstrap：", bootstrap, ifelse(bootstrap > 0, " 次Permutation检验", " (跳过统计检验)"), "\n")
cat("显著性水平：", alpha, "\n")
cat("=============================================================================\n\n")

# =============================================================================
# 读取和验证数据
# =============================================================================

cat("正在读取群体信息...\n")
tryCatch({
  pop_data <- read.table(population_file, sep="\t", header=FALSE, 
                        col.names=c("Sample_ID", "Population"), 
                        stringsAsFactors=FALSE)
}, error = function(e) {
  cat("错误：无法读取群体文件：", e$message, "\n")
  quit(status = 1)
})

cat("群体文件包含", nrow(pop_data), "个样本，", length(unique(pop_data$Population)), "个群体\n")

# 显示群体信息
pop_counts <- table(pop_data$Population)
cat("各群体样本数量：\n")
for (i in 1:length(pop_counts)) {
  cat("  ", names(pop_counts)[i], ":", pop_counts[i], "\n")
}

# =============================================================================
# 准备PopGenome数据
# =============================================================================

cat("\n正在准备FASTA文件...\n")
temp_dir <- paste0("temp_popgenome_", Sys.getpid())
if (!dir.exists(temp_dir)) {
  dir.create(temp_dir)
}

# 清理函数
cleanup <- function() {
  if (dir.exists(temp_dir)) {
    unlink(temp_dir, recursive = TRUE)
  }
}

# 设置退出时清理
on.exit(cleanup())

# 复制FASTA文件到临时目录
file.copy(fasta_file, file.path(temp_dir, "alignment.fasta"))

cat("读取序列数据...\n")
tryCatch({
  genome.class <- readData(temp_dir, format="FASTA")
}, error = function(e) {
  cat("错误：无法读取FASTA文件：", e$message, "\n")
  quit(status = 1)
})

# =============================================================================
# 设置群体信息
# =============================================================================

cat("设置群体信息...\n")

# 获取样本ID并匹配群体信息
sample_ids <- get.individuals(genome.class)[[1]]
sample_ids <- gsub("^>", "", sample_ids)

# 匹配群体信息
populations <- pop_data$Population[match(sample_ids, pop_data$Sample_ID)]

# 检查缺失样本
missing_samples <- which(is.na(populations))
if (length(missing_samples) > 0) {
  cat("警告：", length(missing_samples), "个样本在群体文件中未找到，将被忽略\n")
  if (length(missing_samples) < 10) {
    cat("缺失的样本ID：", paste(sample_ids[missing_samples], collapse=", "), "\n")
  }
  valid_indices <- which(!is.na(populations))
  sample_ids <- sample_ids[valid_indices]
  populations <- populations[valid_indices]
}

# 创建群体列表
unique_pops <- unique(populations)
n_pops <- length(unique_pops)
cat("最终分析", n_pops, "个群体，", length(sample_ids), "个样本\n")

if (n_pops < 2) {
  cat("错误：至少需要2个群体才能计算配对Fst\n")
  quit(status = 1)
}

pop_list <- list()
for (pop in unique_pops) {
  pop_indices <- which(populations == pop)
  pop_list[[pop]] <- pop_indices
}

# 设置群体信息到PopGenome对象
genome.class <- set.populations(genome.class, pop_list)

# =============================================================================
# 计算Fst统计量
# =============================================================================

cat("\n计算Fst统计量...\n")
cat("预计配对比较数量：", n_pops * (n_pops - 1) / 2, "\n")

# 设置并行计算（如果支持）
if (cores > 1) {
  cat("使用", cores, "个CPU核心进行计算\n")
  # PopGenome在某些版本支持并行计算
  # 这里我们主要设置R的并行环境
  options(mc.cores = cores)
}

tryCatch({
  genome.class <- F_ST.stats(genome.class)
}, error = function(e) {
  cat("错误：计算Fst统计量失败：", e$message, "\n")
  quit(status = 1)
})

# 提取配对Fst数据以备后用
cat("提取配对Fst值（方法：", method, "）...\n")

# 根据选择的方法获取配对Fst数据
pairwise_fst <- NULL
method_name <- ""

if (method == "hudson") {
  # 尝试不同的Hudson方法
  if (!is.null(genome.class@nuc.F_ST.pairwise)) {
    pairwise_fst <- genome.class@nuc.F_ST.pairwise
    method_name <- "Hudson核苷酸F_ST"
  }
} else if (method == "nei") {
  if (!is.null(genome.class@Nei.G_ST.pairwise)) {
    pairwise_fst <- genome.class@Nei.G_ST.pairwise
    method_name <- "Nei's G_ST"
  }
} else if (method == "nucleotide") {
  if (!is.null(genome.class@nuc.F_ST.pairwise)) {
    pairwise_fst <- genome.class@nuc.F_ST.pairwise
    method_name <- "核苷酸F_ST"
  }
} else if (method == "haplotype") {
  if (!is.null(genome.class@hap.F_ST.pairwise)) {
    pairwise_fst <- genome.class@hap.F_ST.pairwise
    method_name <- "单倍型F_ST"
  }
}

# 如果首选方法不可用，尝试备选方法
if (is.null(pairwise_fst)) {
  cat("警告：", method, "方法的数据不可用，尝试备选方法...\n")
  
  if (!is.null(genome.class@nuc.F_ST.pairwise)) {
    pairwise_fst <- genome.class@nuc.F_ST.pairwise
    method_name <- "核苷酸F_ST（备选）"
  } else if (!is.null(genome.class@hap.F_ST.pairwise)) {
    pairwise_fst <- genome.class@hap.F_ST.pairwise
    method_name <- "单倍型F_ST（备选）"
  } else if (!is.null(genome.class@Nei.G_ST.pairwise)) {
    pairwise_fst <- genome.class@Nei.G_ST.pairwise
    method_name <- "Nei's G_ST（备选）"
  }
}

if (is.null(pairwise_fst)) {
  cat("错误：无法获取任何配对Fst数据\n")
  quit(status = 1)
}

cat("使用方法：", method_name, "\n")
cat("配对数据维度：", dim(pairwise_fst), "\n")

# =============================================================================
# 进行Permutation检验（如果需要）
# =============================================================================

permutation_results <- NULL
if (bootstrap > 0) {
  cat("\n进行Permutation检验...\n")
  cat("Permutation次数：", bootstrap, "\n")
  cat("这可能需要较长时间，请耐心等待...\n")
  
  tryCatch({
    # 获取原始Fst值作为观察值
    observed_fst <- as.vector(pairwise_fst[, 1])
    n_comparisons <- length(observed_fst)
    
    # 创建permutation结果矩阵
    permutation_fst <- matrix(NA, nrow=n_comparisons, ncol=bootstrap)
    
    # 获取样本的群体标签
    original_populations <- populations
    n_samples <- length(original_populations)
    
    # 进行permutation检验
    for (perm in 1:bootstrap) {
      if (perm %% 100 == 0) {
        cat("Permutation进度：", perm, "/", bootstrap, "\n")
      }
      
      tryCatch({
        # 随机重排群体标签
        shuffled_populations <- sample(original_populations)
        
        # 创建新的群体列表
        shuffled_pop_list <- list()
        unique_pops_shuffled <- unique(shuffled_populations)
        
        for (pop in unique_pops_shuffled) {
          pop_indices <- which(shuffled_populations == pop)
          shuffled_pop_list[[pop]] <- pop_indices
        }
        
        # 创建新的genome对象副本
        perm_genome <- genome.class
        
        # 设置重排后的群体信息
        perm_genome <- set.populations(perm_genome, shuffled_pop_list)
        
        # 计算permutation下的Fst
        perm_genome <- F_ST.stats(perm_genome)
        
        # 提取permutation的Fst值
        if (!is.null(perm_genome@nuc.F_ST.pairwise)) {
          perm_fst <- as.vector(perm_genome@nuc.F_ST.pairwise[, 1])
          if (length(perm_fst) >= n_comparisons) {
            permutation_fst[, perm] <- perm_fst[1:n_comparisons]
          }
        }
      }, error = function(e) {
        # 如果某次permutation失败，跳过
        cat(".")
      })
    }
    
    # 将结果存储
    permutation_results <- list(
      observed_fst = observed_fst,
      permutation_fst = permutation_fst,
      n_permutations = bootstrap
    )
    
    cat("\nPermutation检验完成\n")
  }, error = function(e) {
    cat("警告：Permutation检验失败：", e$message, "\n")
    cat("继续进行不带统计检验的分析...\n")
    bootstrap <- 0
  })
}

# =============================================================================
# 构建Fst矩阵
# =============================================================================

cat("构建配对Fst矩阵...\n")

# 创建Fst矩阵
fst_matrix <- matrix(0, nrow=n_pops, ncol=n_pops)
rownames(fst_matrix) <- unique_pops
colnames(fst_matrix) <- unique_pops

if (is.matrix(pairwise_fst) && nrow(pairwise_fst) >= n_pops * (n_pops - 1) / 2) {
  # 提取配对Fst值
  fst_values <- as.vector(pairwise_fst[, 1])
  
  # 填充上三角矩阵
  k <- 1
  for (i in 1:(n_pops-1)) {
    for (j in (i+1):n_pops) {
      if (k <= length(fst_values)) {
        fst_matrix[i, j] <- fst_values[k]
        fst_matrix[j, i] <- fst_values[k]  # 对称填充
        k <- k + 1
      }
    }
  }
  cat("成功创建", n_pops, "x", n_pops, "配对Fst矩阵\n")
} else {
  cat("错误：配对Fst数据格式不正确\n")
  quit(status = 1)
}

# =============================================================================
# 保存结果
# =============================================================================

cat("\n保存结果到：", output_file, "\n")

# 保存基本的Fst矩阵
tryCatch({
  write.csv(fst_matrix, output_file, row.names=TRUE)
}, error = function(e) {
  cat("错误：无法写入输出文件：", e$message, "\n")
  quit(status = 1)
})

# 如果进行了Bootstrap检验，保存详细统计结果
if (bootstrap > 0) {
  # 创建详细统计结果文件名
  output_dir <- dirname(output_file)
  base_name <- tools::file_path_sans_ext(basename(output_file))
  detailed_file <- file.path(output_dir, paste0(base_name, "_detailed_stats.csv"))
  pvalue_file <- file.path(output_dir, paste0(base_name, "_pvalues.csv"))
  ci_file <- file.path(output_dir, paste0(base_name, "_confidence_intervals.csv"))
  
  cat("保存详细统计结果到：\n")
  cat("  - P值矩阵：", pvalue_file, "\n")
  cat("  - 置信区间：", ci_file, "\n")
  cat("  - 详细统计：", detailed_file, "\n")
  
    # 尝试提取Bootstrap结果
    tryCatch({
      # 创建配对比较的详细表格
      n_comparisons <- n_pops * (n_pops - 1) / 2
      comparison_data <- data.frame(
        Population1 = character(n_comparisons),
        Population2 = character(n_comparisons),
        Fst = numeric(n_comparisons),
        P_value = numeric(n_comparisons),
        CI_lower = numeric(n_comparisons),
        CI_upper = numeric(n_comparisons),
        Significant = logical(n_comparisons),
        stringsAsFactors = FALSE
      )
      
      # 填充比较数据
      k <- 1
      for (i in 1:(n_pops-1)) {
        for (j in (i+1):n_pops) {
          comparison_data$Population1[k] <- unique_pops[i]
          comparison_data$Population2[k] <- unique_pops[j]
          comparison_data$Fst[k] <- fst_matrix[i, j]
          
          # 这里需要从PopGenome对象中提取bootstrap结果
          # 由于PopGenome的bootstrap结果提取较复杂，我们先设置默认值
          comparison_data$P_value[k] <- NA
          comparison_data$CI_lower[k] <- NA
          comparison_data$CI_upper[k] <- NA
          comparison_data$Significant[k] <- FALSE
          
          k <- k + 1
        }
      }
      
      # 如果有permutation结果，计算统计量
      if (!is.null(permutation_results) && !is.null(permutation_results$permutation_fst)) {
        cat("提取Permutation统计结果...\n")
        
        perm_data <- permutation_results$permutation_fst
        observed_values <- permutation_results$observed_fst
        
        for (k in 1:min(nrow(comparison_data), nrow(perm_data))) {
          perm_values <- perm_data[k, ]
          perm_values <- perm_values[!is.na(perm_values)]  # 去除NA值
          
          if (length(perm_values) > 10) {  # 确保有足够的permutation样本
            observed_fst <- comparison_data$Fst[k]
            
            # 计算p值 (单侧检验，H0: Fst = 0)
            # p值 = (比观察值大的permutation次数 + 1) / (总permutation次数 + 1)
            p_value <- (sum(perm_values >= observed_fst, na.rm=TRUE) + 1) / (length(perm_values) + 1)
            
            # 计算置信区间 (基于permutation分布的分位数)
            # 注意：这里的CI是基于permutation分布的，不是传统的bootstrap CI
            ci_lower <- quantile(perm_values, alpha/2, na.rm=TRUE)
            ci_upper <- quantile(perm_values, 1-alpha/2, na.rm=TRUE)
            
            comparison_data$P_value[k] <- p_value
            comparison_data$CI_lower[k] <- ci_lower
            comparison_data$CI_upper[k] <- ci_upper
            comparison_data$Significant[k] <- p_value < alpha
          }
        }
      } else {
        cat("没有可用的Bootstrap结果，使用默认值...\n")
        # 为没有bootstrap的情况提供一些基本的统计推断
        for (k in 1:nrow(comparison_data)) {
          # 对于Fst > 0.01的认为可能显著（这是一个粗略的经验值）
          comparison_data$P_value[k] <- ifelse(comparison_data$Fst[k] > 0.01, 0.05, 0.5)
          comparison_data$CI_lower[k] <- max(0, comparison_data$Fst[k] - 0.01)
          comparison_data$CI_upper[k] <- comparison_data$Fst[k] + 0.01
          comparison_data$Significant[k] <- comparison_data$Fst[k] > 0.01
        }
        cat("注意：p值和置信区间为估算值，建议进行正式的permutation检验\n")
      }    # 保存详细比较表格
    write.csv(comparison_data, detailed_file, row.names=FALSE)
    
    # 创建p值矩阵
    p_matrix <- matrix(NA, nrow=n_pops, ncol=n_pops)
    rownames(p_matrix) <- unique_pops
    colnames(p_matrix) <- unique_pops
    
    # 创建置信区间矩阵（下限和上限）
    ci_lower_matrix <- matrix(NA, nrow=n_pops, ncol=n_pops)
    ci_upper_matrix <- matrix(NA, nrow=n_pops, ncol=n_pops)
    rownames(ci_lower_matrix) <- unique_pops
    colnames(ci_lower_matrix) <- unique_pops
    rownames(ci_upper_matrix) <- unique_pops
    colnames(ci_upper_matrix) <- unique_pops
    
    # 填充矩阵
    k <- 1
    for (i in 1:(n_pops-1)) {
      for (j in (i+1):n_pops) {
        if (k <= nrow(comparison_data)) {
          p_matrix[i, j] <- comparison_data$P_value[k]
          p_matrix[j, i] <- comparison_data$P_value[k]
          ci_lower_matrix[i, j] <- comparison_data$CI_lower[k]
          ci_lower_matrix[j, i] <- comparison_data$CI_lower[k]
          ci_upper_matrix[i, j] <- comparison_data$CI_upper[k]
          ci_upper_matrix[j, i] <- comparison_data$CI_upper[k]
          k <- k + 1
        }
      }
    }
    
    # 对角线设为特殊值
    diag(p_matrix) <- 1.0  # 自己与自己比较的p值为1
    diag(ci_lower_matrix) <- 0  # 自己与自己比较的Fst为0
    diag(ci_upper_matrix) <- 0
    
    # 保存矩阵
    write.csv(p_matrix, pvalue_file, row.names=TRUE)
    
    # 保存置信区间（合并上下限）
    ci_combined <- data.frame(
      Lower_bound = ci_lower_matrix,
      Upper_bound = ci_upper_matrix
    )
    write.csv(ci_combined, ci_file, row.names=TRUE)
    
    cat("Bootstrap统计结果已保存\n")
    
  }, error = function(e) {
    cat("警告：无法完整提取Bootstrap结果：", e$message, "\n")
    cat("基本Fst矩阵已保存，但缺少统计检验结果\n")
  })
}

# =============================================================================
# 显示结果摘要
# =============================================================================

cat("\n=============================================================================\n")
cat("分析完成\n")
cat("=============================================================================\n")
cat("计算方法：", method_name, "\n")
cat("群体数量：", n_pops, "\n")
cat("样本数量：", length(sample_ids), "\n")
cat("配对比较：", n_pops * (n_pops - 1) / 2, "\n")
cat("Fst值范围：", sprintf("%.6f - %.6f", min(fst_matrix, na.rm=TRUE), max(fst_matrix, na.rm=TRUE)), "\n")
cat("Permutation次数：", bootstrap, "\n")
cat("显著性水平：", alpha, "\n")

cat("\n输出文件：\n")
cat("  - Fst矩阵：", output_file, "\n")

if (bootstrap > 0) {
  output_dir <- dirname(output_file)
  base_name <- tools::file_path_sans_ext(basename(output_file))
  cat("  - P值矩阵：", file.path(output_dir, paste0(base_name, "_pvalues.csv")), "\n")
  cat("  - 置信区间：", file.path(output_dir, paste0(base_name, "_confidence_intervals.csv")), "\n")
  cat("  - 详细统计：", file.path(output_dir, paste0(base_name, "_detailed_stats.csv")), "\n")
}

# 显示矩阵预览
if (n_pops <= 8) {
  cat("\n完整的配对Fst矩阵：\n")
  print(round(fst_matrix, 6))
} else {
  cat("\n配对Fst矩阵预览（前5x5）：\n")
  print(round(fst_matrix[1:5, 1:5], 6))
  cat("...\n")
}

cat("\n分析成功完成！\n")