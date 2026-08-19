############################################################
# GSE185062 miRNA analysis
#
# 主分析：
# Control vs NAFL
# 目标：
# 1. DESeq2 差异 miRNA
# 2. 先筛选高表达量且高AUC的miRNA，再搜索最佳三miRNA panel
# 3. Logistic panel + UDR-like score
# 4. ROC, CV ROC, 混淆矩阵, 热图, 箱线图
#
# 补充分析：
# 1. F0_NAFLD 内部真实 NAS 相关性
# 2. Control + F0_NAFLD NAS 分组展示
#
# 特别设置：
# miR-3687 完全排除
#
# 本修订版分析顺序：
# 1. 先删除CTRL_/HK_技术对照、内参及人工排除的miRNA
# 2. 评价真实miRNA的原始read count丰度、检出率和组内稳定性
# 3. 先用“高表达量 + 高单miRNA AUC”作为硬性预筛条件
# 4. DESeq2的padj和log2FC保留为生物学支持证据，不再作为首轮硬筛条件
# 5. 从预筛后的候选池中遍历3-miRNA组合，并用重复分层交叉验证评价
# 6. CV-AUC近似时，优先选择表达量和检出率更高的组合
############################################################

suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pROC)
  library(ggrepel)
  library(pheatmap)
})

set.seed(123)
select <- dplyr::select
filter <- dplyr::filter
arrange <- dplyr::arrange
mutate <- dplyr::mutate
summarise <- dplyr::summarise
left_join <- dplyr::left_join
pull <- dplyr::pull
############################################################
# 0. 参数区
############################################################


out_dir <- "C:/Users/24791/Desktop/miRNA/GEO-NAFLD/output"

counts_file <- file.path(out_dir, "GSE185062_miRNA_counts_matrix.csv")
meta_file   <- file.path(out_dir, "GSE185062_sample_info.csv")

final_dir <- file.path(
  out_dir,
  "FINAL_high_expression_high_AUC_then_panel_no_controls"
)

dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

fig_dir <- file.path(final_dir, "figures")
tab_dir <- file.path(final_dir, "tables")
nas_dir <- file.path(final_dir, "NAS_analysis")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(nas_dir, recursive = TRUE, showWarnings = FALSE)

min_total_count <- 300
min_base_mean   <- 10

padj_cut <- 0.05
lfc_cut  <- 1

top_heatmap_n <- 30

k_folds   <- 5
n_repeats <- 30

exclude_miRNAs <- c("miR-3687")

# 技术对照、外源对照和管家/内参基因不参与miRNA biomarker筛选
# 保留到单独QC表，但在PCA、DESeq2、AUC和panel建模前移除
exclude_feature_regex <- "^(CTRL_|HK_)"


write_csv(
  tibble(excluded_miRNA = exclude_miRNAs),
  file.path(tab_dir, "excluded_miRNAs.csv")
)
############################################################
# miRNA原始表达丰度与高AUC预筛参数
############################################################

# 原始read count达到该值，视为该样本中可检出
min_detect_count <- 10

# expression_score内部权重：原始丰度权重最高
expression_abundance_weight <- 0.50
expression_detection_weight <- 0.35
expression_stability_weight <- 0.15

# 第一阶段：高表达量硬筛条件
# 这些阈值评价检测可行性，而不是用于DESeq2归一化
high_expr_min_mean_raw_count <- 100
high_expr_min_median_raw_count <- 20
high_expr_min_detection_rate <- 0.80

# 为避免某一个组几乎完全无法检出，要求两组中较低检出率至少达到该值
# 下调型marker在疾病组中可能较低，因此这里不宜设得过高
high_expr_min_group_detection_rate <- 0.50

# 第二阶段：单miRNA诊断能力硬筛条件
high_auc_cut <- 0.80

# 通过硬筛后，候选排序只综合表达可检出性和单miRNA AUC
preselect_expression_weight <- 0.50
preselect_auc_weight <- 0.50

stopifnot(
  abs(
    expression_abundance_weight +
      expression_detection_weight +
      expression_stability_weight - 1
  ) < 1e-8,
  abs(
    preselect_expression_weight +
      preselect_auc_weight - 1
  ) < 1e-8
)

############################################################
# 双功能panel筛选参数
############################################################

# Control vs NAFL
dual_disease_padj_cut <- 0.05
dual_disease_lfc_cut  <- 1
dual_disease_auc_cut  <- 0.75

# NAS 1-2 vs NAS 3-4
dual_nas_padj_cut <- 0.05
dual_nas_lfc_cut  <- 0.5
dual_nas_auc_cut  <- 0.70

# panel中任意两个miRNA的最大相关性
dual_cor_cut <- 0.85

# 最多使用排名前多少个共同候选进行三组合搜索
# 12个候选共220种三组合，计算量相对可控
dual_candidate_n <- 12

# 三组合初筛时重复交叉验证次数
dual_cv_repeats <- 10


############################################################
# 疾病诊断3-miRNA panel筛选参数
############################################################

# DESeq2支持证据阈值：仅用于结果标记和辅助解释，不作为首轮硬筛条件
panel_padj_cut <- 0.05
panel_lfc_cut  <- 1

# panel候选的单miRNA AUC硬阈值由high_auc_cut控制
panel_auc_cut <- high_auc_cut

# panel内任意两个miRNA允许的最大绝对Spearman相关系数
panel_cor_cut <- 0.85

# 最多选择综合排名前多少个miRNA进行三组合搜索
# 12个miRNA一共有220种三组合
panel_candidate_n <- 12

# 三组合筛选使用的重复交叉验证次数
panel_cv_repeats <- 30
############################################################
# 1. 工具函数
############################################################

save_gg <- function(p, filename, width = 8, height = 6) {
  ggsave(
    file.path(fig_dir, paste0(filename, ".png")),
    p,
    width = width,
    height = height,
    dpi = 300
  )
  ggsave(
    file.path(fig_dir, paste0(filename, ".pdf")),
    p,
    width = width,
    height = height
  )
}

make_stratified_folds <- function(y, k = 5) {
  y <- as.factor(y)
  folds <- vector("list", k)
  
  for (cls in levels(y)) {
    idx <- which(y == cls)
    idx <- sample(idx)
    fold_id <- rep(seq_len(k), length.out = length(idx))
    
    for (i in seq_len(k)) {
      folds[[i]] <- c(folds[[i]], idx[fold_id == i])
    }
  }
  
  folds
}

scale_train_test <- function(train_df, test_df, vars) {
  for (v in vars) {
    m <- mean(train_df[[v]], na.rm = TRUE)
    s <- sd(train_df[[v]], na.rm = TRUE)
    
    if (is.na(s) || s == 0) s <- 1
    
    train_df[[v]] <- (train_df[[v]] - m) / s
    test_df[[v]]  <- (test_df[[v]]  - m) / s
  }
  
  list(train = train_df, test = test_df)
}

run_deseq2 <- function(count_mat, col_df, design_formula, contrast_vec) {
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData   = col_df,
    design    = design_formula
  )
  
  dds <- dds[rowSums(counts(dds)) >= min_total_count, ]
  
  dds <- DESeq(dds)
  
  res <- results(dds, contrast = contrast_vec)
  
  res_df <- as.data.frame(res) %>%
    rownames_to_column("miRNA") %>%
    arrange(padj)
  
  list(dds = dds, res = res_df)
}

get_log_norm <- function(count_mat, col_df) {
  dds <- DESeqDataSetFromMatrix(
    countData = round(count_mat),
    colData   = col_df,
    design    = ~ 1
  )
  
  dds <- dds[rowSums(counts(dds)) >= min_total_count, ]
  
  dds <- estimateSizeFactors(dds)
  
  norm_counts <- counts(dds, normalized = TRUE)
  
  log2(norm_counts + 1)
}

select_2up_1down <- function(res_df, exclude = exclude_miRNAs) {
  candidate <- res_df %>%
    filter(!miRNA %in% exclude) %>%
    filter(!is.na(pvalue)) %>%
    mutate(
      p_rank = if_else(is.na(pvalue), 1, pvalue),
      rank_score = -log10(p_rank + 1e-300) + log10(baseMean + 1)
    )
  
  up <- candidate %>%
    filter(
      log2FoldChange > 0,
      baseMean >= min_base_mean
    ) %>%
    arrange(desc(rank_score))
  
  down <- candidate %>%
    filter(
      log2FoldChange < 0,
      baseMean >= min_base_mean
    ) %>%
    arrange(desc(rank_score))
  
  if (nrow(up) < 2) {
    stop("上调候选 miRNA 少于 2 个，请降低 min_base_mean 或检查分组。")
  }
  
  if (nrow(down) < 1) {
    stop("下调候选 miRNA 少于 1 个，请降低 min_base_mean 或检查分组。")
  }
  
  selected <- c(
    up %>% slice_head(n = 2) %>% pull(miRNA),
    down %>% slice_head(n = 1) %>% pull(miRNA)
  )
  
  list(
    selected = selected,
    up = up,
    down = down,
    candidate = candidate
  )
}

roc_table_one <- function(response, predictor, marker_name) {
  roc_obj <- roc(response, predictor, quiet = TRUE)
  
  tibble(
    marker = marker_name,
    AUC = as.numeric(auc(roc_obj)),
    CI_low = as.numeric(ci.auc(roc_obj)[1]),
    CI_mid = as.numeric(ci.auc(roc_obj)[2]),
    CI_high = as.numeric(ci.auc(roc_obj)[3])
  )
}

roc_curve_df <- function(response, predictor, marker_name) {
  roc_obj <- roc(response, predictor, quiet = TRUE)
  
  tibble(
    marker = marker_name,
    specificity = roc_obj$specificities,
    sensitivity = roc_obj$sensitivities
  )
}

make_confusion_from_roc <- function(response, predictor, label = "panel") {
  roc_obj <- roc(response, predictor, quiet = TRUE)
  
  best <- coords(
    roc_obj,
    x = "best",
    best.method = "youden",
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  threshold <- as.numeric(best$threshold)
  
  pred_class <- if_else(predictor >= threshold, 1, 0)
  
  tab <- table(
    truth = factor(response, levels = c(0, 1)),
    pred  = factor(pred_class, levels = c(0, 1))
  )
  
  tibble(
    model = label,
    threshold = threshold,
    TN = tab["0", "0"],
    FP = tab["0", "1"],
    FN = tab["1", "0"],
    TP = tab["1", "1"],
    sensitivity = as.numeric(best$sensitivity),
    specificity = as.numeric(best$specificity),
    accuracy = (tab["0", "0"] + tab["1", "1"]) / sum(tab)
  )
}
############################################################
# miRNA原始表达量、检出率及组内稳定性评价
############################################################

# 安全的百分位排名：最小值接近0，最大值为1
safe_percent_rank <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  n_ok <- sum(ok)

  if (n_ok == 0) {
    return(out)
  }

  if (n_ok == 1) {
    out[ok] <- 1
    return(out)
  }

  out[ok] <- (
    rank(x[ok], ties.method = "average") - 1
  ) / (n_ok - 1)

  out
}

# 安全计算变异系数；均值为0时返回NA
safe_cv <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]

  if (length(x) < 2) {
    return(NA_real_)
  }

  current_mean <- mean(x)

  if (!is.finite(current_mean) || current_mean <= 0) {
    return(NA_real_)
  }

  current_sd <- sd(x)

  if (!is.finite(current_sd)) {
    return(NA_real_)
  }

  current_sd / current_mean
}

calculate_expression_score <- function(
    count_matrix,
    meta_df = NULL,
    group_col = "group",
    detect_count_cutoff = min_detect_count
) {

  count_matrix <- as.matrix(count_matrix)
  storage.mode(count_matrix) <- "numeric"

  if (is.null(rownames(count_matrix))) {
    stop("count_matrix必须包含miRNA行名。")
  }

  if (is.null(colnames(count_matrix))) {
    stop("count_matrix必须包含样本列名。")
  }

  # 全部样本的原始表达评价
  mean_raw_count <- rowMeans(
    count_matrix,
    na.rm = TRUE
  )

  median_raw_count <- apply(
    count_matrix,
    1,
    median,
    na.rm = TRUE
  )

  detection_rate <- rowMeans(
    count_matrix >= detect_count_cutoff,
    na.rm = TRUE
  )

  overall_expression_cv <- apply(
    count_matrix,
    1,
    safe_cv
  )

  # 默认先建立空的分组指标
  mean_raw_count_Control <- rep(NA_real_, nrow(count_matrix))
  mean_raw_count_NAFL <- rep(NA_real_, nrow(count_matrix))
  median_raw_count_Control <- rep(NA_real_, nrow(count_matrix))
  median_raw_count_NAFL <- rep(NA_real_, nrow(count_matrix))
  detection_rate_Control <- rep(NA_real_, nrow(count_matrix))
  detection_rate_NAFL <- rep(NA_real_, nrow(count_matrix))
  cv_Control <- rep(NA_real_, nrow(count_matrix))
  cv_NAFL <- rep(NA_real_, nrow(count_matrix))

  # 若提供meta，则分别计算Control和NAFL内部指标
  if (!is.null(meta_df)) {
    if (is.null(rownames(meta_df))) {
      stop("meta_df必须使用样本名作为行名。")
    }

    if (!group_col %in% colnames(meta_df)) {
      stop(
        paste0(
          "meta_df中不存在分组列：",
          group_col
        )
      )
    }

    common_samples <- intersect(
      colnames(count_matrix),
      rownames(meta_df)
    )

    if (length(common_samples) == 0) {
      stop("count_matrix与meta_df之间没有共同样本。")
    }

    meta_use <- meta_df[
      common_samples,
      ,
      drop = FALSE
    ]

    group_value <- as.character(
      meta_use[[group_col]]
    )

    control_samples <- rownames(meta_use)[
      group_value == "Control"
    ]

    nafl_samples <- rownames(meta_use)[
      group_value == "NAFL"
    ]

    if (length(control_samples) > 0) {
      control_mat <- count_matrix[
        ,
        control_samples,
        drop = FALSE
      ]

      mean_raw_count_Control <- rowMeans(
        control_mat,
        na.rm = TRUE
      )

      median_raw_count_Control <- apply(
        control_mat,
        1,
        median,
        na.rm = TRUE
      )

      detection_rate_Control <- rowMeans(
        control_mat >= detect_count_cutoff,
        na.rm = TRUE
      )

      cv_Control <- apply(
        control_mat,
        1,
        safe_cv
      )
    }

    if (length(nafl_samples) > 0) {
      nafl_mat <- count_matrix[
        ,
        nafl_samples,
        drop = FALSE
      ]

      mean_raw_count_NAFL <- rowMeans(
        nafl_mat,
        na.rm = TRUE
      )

      median_raw_count_NAFL <- apply(
        nafl_mat,
        1,
        median,
        na.rm = TRUE
      )

      detection_rate_NAFL <- rowMeans(
        nafl_mat >= detect_count_cutoff,
        na.rm = TRUE
      )

      cv_NAFL <- apply(
        nafl_mat,
        1,
        safe_cv
      )
    }
  }

  # 使用组内CV评价稳定性，避免真实组间差异抬高总体CV
  within_group_cv_matrix <- cbind(
    cv_Control,
    cv_NAFL
  )

  mean_within_group_cv <- apply(
    within_group_cv_matrix,
    1,
    function(x) {
      x <- x[is.finite(x)]
      if (length(x) == 0) NA_real_ else mean(x)
    }
  )

  # 两组中较低的检出率只作为结果输出，不直接作为硬性过滤
  min_group_detection_rate <- apply(
    cbind(
      detection_rate_Control,
      detection_rate_NAFL
    ),
    1,
    function(x) {
      x <- x[is.finite(x)]
      if (length(x) == 0) NA_real_ else min(x)
    }
  )

  expr_summary <- tibble(
    miRNA = rownames(count_matrix),
    mean_raw_count = mean_raw_count,
    median_raw_count = median_raw_count,
    detection_rate = detection_rate,
    overall_expression_cv = overall_expression_cv,
    mean_raw_count_Control = mean_raw_count_Control,
    mean_raw_count_NAFL = mean_raw_count_NAFL,
    median_raw_count_Control = median_raw_count_Control,
    median_raw_count_NAFL = median_raw_count_NAFL,
    detection_rate_Control = detection_rate_Control,
    detection_rate_NAFL = detection_rate_NAFL,
    min_group_detection_rate = min_group_detection_rate,
    cv_Control = cv_Control,
    cv_NAFL = cv_NAFL,
    expression_cv = mean_within_group_cv
  ) %>%
    mutate(
      # 原始表达越高，丰度排名越高
      abundance_rank = safe_percent_rank(
        log10(mean_raw_count + 1)
      ),

      # 在更多样本中达到检测阈值，检出排名越高
      detection_rank = safe_percent_rank(
        detection_rate
      ),

      # 组内CV越小，稳定性排名越高
      stability_rank = safe_percent_rank(
        -expression_cv
      ),

      # 表达可检出性综合评分：丰度权重最高
      expression_score =
        expression_abundance_weight * abundance_rank +
        expression_detection_weight * detection_rank +
        expression_stability_weight * stability_rank
    ) %>%
    arrange(
      desc(expression_score),
      desc(mean_raw_count),
      desc(detection_rate)
    )

  expr_summary
}

# 汇总一个panel的表达可检出性，用于近似最佳panel之间的优先选择
summarize_panel_expression <- function(
    features,
    expression_df
) {

  current_expression <- expression_df %>%
    filter(
      miRNA %in% features
    )

  if (nrow(current_expression) != length(features)) {
    missing_features <- setdiff(
      features,
      current_expression$miRNA
    )

    stop(
      paste0(
        "表达评价表中缺少以下miRNA：",
        paste(missing_features, collapse = ", ")
      )
    )
  }

  tibble(
    panel_mean_expression_score = mean(
      current_expression$expression_score,
      na.rm = TRUE
    ),
    panel_min_expression_score = min(
      current_expression$expression_score,
      na.rm = TRUE
    ),
    panel_mean_raw_count = mean(
      current_expression$mean_raw_count,
      na.rm = TRUE
    ),
    panel_median_raw_count = median(
      current_expression$median_raw_count,
      na.rm = TRUE
    ),
    panel_mean_detection_rate = mean(
      current_expression$detection_rate,
      na.rm = TRUE
    ),
    panel_min_detection_rate = min(
      current_expression$detection_rate,
      na.rm = TRUE
    )
  )
}

############################################################
# 2. 读取数据
############################################################

counts_raw <- read_csv(counts_file, show_col_types = FALSE) %>%
  column_to_rownames("Gene")

meta <- read_csv(meta_file, show_col_types = FALSE)

meta <- meta %>%
  filter(GSM %in% colnames(counts_raw))

counts_raw <- counts_raw[, meta$GSM]

meta <- meta %>%
  column_to_rownames("GSM")

meta <- meta[colnames(counts_raw), , drop = FALSE]

stopifnot(all(colnames(counts_raw) == rownames(meta)))

meta <- meta %>%
  mutate(
    fibrosis_stage = suppressWarnings(as.numeric(fibrosis_stage)),
    nas            = suppressWarnings(as.numeric(nas)),
    saf_activity   = suppressWarnings(as.numeric(saf_activity))
  )
############################################################
# 仅保留 group = Control 和 NAFL
############################################################

meta <- meta %>%
  filter(
    group %in% c("Control", "NAFL"),
    !(group == "NAFL" & nas == 5)
  )

counts_raw <- counts_raw[, rownames(meta), drop = FALSE]


stopifnot(all(colnames(counts_raw) == rownames(meta)))


meta <- meta %>%
  mutate(
    group = factor(
      group,
      levels = c("Control", "NAFL")
    ),
    outcome = if_else(group == "NAFL", 1, 0)
  )


cat("\n保留 Control 和 NAFL 后样本量：\n")
print(table(meta$group))

cat("\n原始样本分组：\n")
print(table(meta$case_control, useNA = "ifany"))

cat("\nNAFLD fibrosis_stage 分布：\n")
print(table(meta$fibrosis_stage[meta$case_control == "NAFLD"], useNA = "ifany"))

cat("\nF0_NAFLD 真实 NAS 分布：\n")
print(table(meta$nas[meta$case_control == "NAFLD" & meta$fibrosis_stage == 0], useNA = "ifany"))

############################################################
# 3. 保留原始 miRNA（不合并3p/5p）
############################################################

counts_total <- round(counts_raw)

# 识别技术对照、管家基因以及人工排除的miRNA
feature_names <- rownames(counts_total)

technical_control_flag <- stringr::str_detect(
  feature_names,
  stringr::regex(
    exclude_feature_regex,
    ignore_case = TRUE
  )
)

manual_exclusion_flag <- feature_names %in% exclude_miRNAs

excluded_feature_flag <- technical_control_flag | manual_exclusion_flag

excluded_feature_counts <- counts_total[
  excluded_feature_flag,
  ,
  drop = FALSE
]

write_csv(
  as.data.frame(excluded_feature_counts) %>%
    rownames_to_column("feature"),
  file.path(
    tab_dir,
    "Technical_HK_and_manually_excluded_features_QC_only.csv"
  )
)

write_csv(
  tibble(
    feature = feature_names[excluded_feature_flag],
    exclusion_reason = case_when(
      technical_control_flag[excluded_feature_flag] ~ "CTRL_or_HK_feature",
      manual_exclusion_flag[excluded_feature_flag] ~ "manual_exclusion",
      TRUE ~ "other"
    )
  ),
  file.path(
    tab_dir,
    "Excluded_features_before_biomarker_analysis.csv"
  )
)

# 后续PCA、DESeq2、AUC和panel仅使用真实候选miRNA
counts_total <- counts_total[
  !excluded_feature_flag,
  ,
  drop = FALSE
]

cat(
  "\n移除CTRL_/HK_和人工排除项后，进入生物标志物分析的miRNA数量：",
  nrow(counts_total),
  "\n"
)

write_csv(
  as.data.frame(counts_total) %>%
    rownames_to_column("miRNA"),
  file.path(
    tab_dir,
    "raw_biological_miRNA_counts_after_control_removal.csv"
  )
)
############################################################
#PCA检测并剔除异常样本
#
# 目的：
# 在差异分析前检测样本整体表达异常
# 避免异常样本影响DESeq2结果
############################################################


# 构建临时DESeq2对象用于PCA

dds_pca <- DESeqDataSetFromMatrix(
  countData = counts_total,
  colData = meta,
  design = ~ group
)


# 去除低表达miRNA

dds_pca <- dds_pca[
  rowSums(counts(dds_pca)) >= min_total_count,
]


# VST标准化

vsd_pca <- vst(
  dds_pca,
  blind = TRUE
)


# PCA分析

pca_data <- plotPCA(
  vsd_pca,
  intgroup = "group",
  returnData = TRUE
)


percentVar <- round(
  100 * attr(
    pca_data,
    "percentVar"
  )
)



# PCA可视化

p_pca_before <- ggplot(
  pca_data,
  aes(
    PC1,
    PC2,
    color = group,
    label = name
  )
)+
  geom_point(
    size = 3
  )+
  geom_text(
    size = 3,
    vjust=-0.5
  )+
  labs(
    title="PCA before outlier removal",
    x=paste0("PC1: ",percentVar[1],"%"),
    y=paste0("PC2: ",percentVar[2],"%")
  )+
  theme_bw()


print(p_pca_before)


ggsave(
  file.path(
    fig_dir,
    "PCA_before_outlier_removal.png"
  ),
  p_pca_before,
  width = 7,
  height = 6,
  dpi = 300
)



############################################################
# PCA距离计算
############################################################


pca_matrix <- pca_data %>%
  select(
    PC1,
    PC2
  )


pca_center <- colMeans(
  pca_matrix
)


pca_distance <- sqrt(
  rowSums(
    sweep(
      pca_matrix,
      2,
      pca_center
    )^2
  )
)


pca_outlier <- tibble(
  GSM = pca_data$name,
  PC1 = pca_data$PC1,
  PC2 = pca_data$PC2,
  distance = pca_distance
)



# 使用95%分位作为异常阈值

pca_cutoff <- quantile(
  pca_outlier$distance,
  0.95
)


remove_samples <- pca_outlier %>%
  filter(
    distance > pca_cutoff
  ) %>%
  pull(GSM)



cat("\nPCA检测异常样本:\n")
print(remove_samples)



# 保存异常检测结果

write_csv(
  pca_outlier,
  file.path(
    tab_dir,
    "PCA_outlier_detection_results.csv"
  )
)


write_csv(
  tibble(
    removed_GSM = remove_samples
  ),
  file.path(
    tab_dir,
    "PCA_removed_samples.csv"
  )
)



############################################################
# 剔除异常样本
############################################################


if(length(remove_samples)>0){
  
  meta <- meta[
    !rownames(meta) %in% remove_samples,
    ,
    drop = FALSE
  ]
  
  
  counts_total <- counts_total[
    ,
    colnames(counts_total) %in% rownames(meta),
    drop = FALSE
  ]
  
}


stopifnot(
  all(
    colnames(counts_total)==rownames(meta)
  )
)


cat("\nPCA剔除异常后样本量:\n")

print(
  table(meta$group)
)
############################################################
# 4. 构建 Control vs NAFL 主分析数据
############################################################


meta_f0 <- meta %>%
  rownames_to_column("GSM") %>%
  filter(group %in% c("Control", "NAFL")) %>%
  mutate(
    group = factor(
      group,
      levels = c("Control", "NAFL")
    ),
    outcome = if_else(group == "NAFL", 1, 0)
  ) %>%
  column_to_rownames("GSM")


counts_f0 <- counts_total[
  ,
  rownames(meta_f0),
  drop = FALSE
]


stopifnot(
  all(colnames(counts_f0) == rownames(meta_f0))
)


cat("\nControl vs NAFL 样本量：\n")
print(table(meta_f0$group))


write_csv(
  meta_f0 %>% rownames_to_column("GSM"),
  file.path(tab_dir,
            "sample_info_Control_vs_NAFL.csv")
)
############################################################
# 5. DESeq2 差异分析
############################################################

de_f0 <- run_deseq2(
  count_mat = counts_f0,
  col_df = meta_f0,
  design_formula = ~ group,
  contrast_vec = c("group", "NAFL", "Control")
)

dds_f0 <- de_f0$dds
res_f0_df <- de_f0$res

write_csv(
  res_f0_df,
  file.path(tab_dir, "DESeq2_Control_vs_NAFL_results.csv")
)

expression_summary <- calculate_expression_score(
  count_matrix = counts_f0,
  meta_df = meta_f0,
  group_col = "group",
  detect_count_cutoff = min_detect_count
)

required_expression_columns <- c(
  "miRNA",
  "mean_raw_count",
  "median_raw_count",
  "detection_rate",
  "expression_cv",
  "expression_score"
)

if (!all(required_expression_columns %in% colnames(expression_summary))) {
  stop(
    paste0(
      "expression_summary缺少必要列：",
      paste(
        setdiff(
          required_expression_columns,
          colnames(expression_summary)
        ),
        collapse = ", "
      )
    )
  )
}

write_csv(
  expression_summary,
  file.path(
    tab_dir,
    "miRNA_expression_abundance_evaluation.csv"
  )
)

cat("\n原始表达丰度与可检出性评分前20：\n")
print(
  expression_summary %>%
    slice_head(n = 20),
  n = 20,
  width = Inf
)
cat("\nControl vs NAFL 差异分析前 20：\n")
print(head(res_f0_df, 20))
############################################################
# 生成log2标准化表达矩阵
# 用于箱线图展示
############################################################


log_expr_f0 <- get_log_norm(
  counts_f0,
  meta_f0
)


# 去除排除miRNA

log_expr_f0 <- log_expr_f0[
  !rownames(log_expr_f0) %in% exclude_miRNAs,
  ,
  drop = FALSE
]


cat(
  "\n标准化表达矩阵维度:\n"
)

print(
  dim(log_expr_f0)
)
############################################################
# Top20差异miRNA表达箱线图
# Control vs NAFL
############################################################


# 取padj最小的前20个miRNA

top20_miRNA <- res_f0_df %>%
  filter(!is.na(padj)) %>%
  filter(!miRNA %in% exclude_miRNAs) %>%
  arrange(padj) %>%
  slice_head(n = 20) %>%
  pull(miRNA)


cat("\nTop20 miRNA:\n")
print(top20_miRNA)



# 提取表达矩阵

box_expr <- log_expr_f0[
  top20_miRNA,
  ,
  drop = FALSE
]



# 转成长格式

box_df <- as.data.frame(box_expr) %>%
  rownames_to_column("miRNA") %>%
  pivot_longer(
    cols = -miRNA,
    names_to = "GSM",
    values_to = "expression"
  ) %>%
  left_join(
    meta_f0 %>%
      rownames_to_column("GSM") %>%
      select(
        GSM,
        group
      ),
    by="GSM"
  )



# 固定顺序
box_df$miRNA <- factor(
  box_df$miRNA,
  levels = top20_miRNA
)



############################################################
# 绘制箱线图
############################################################


p_box20 <- ggplot(
  box_df,
  aes(
    x = group,
    y = expression,
    fill = group
  )
)+
  
  geom_boxplot(
    width=0.6,
    outlier.shape = NA,
    alpha=0.75
  )+
  
  geom_jitter(
    width=0.15,
    size=1.5,
    alpha=0.7
  )+
  
  facet_wrap(
    ~miRNA,
    scales="free_y",
    ncol=5
  )+
  
  labs(
    title="Top 20 Differentially Expressed miRNAs",
    subtitle="NAFL vs Control",
    x=NULL,
    y="log2 normalized expression",
    fill="Group"
  )+
  
  theme_bw(base_size=13)+
  
  theme(
    plot.title=element_text(
      hjust=0.5,
      face="bold"
    ),
    plot.subtitle=element_text(
      hjust=0.5
    ),
    strip.text=element_text(
      face="bold"
    ),
    legend.position="right",
    panel.grid.minor=element_blank()
  )



print(p_box20)



save_gg(
  p_box20,
  "Top20_DE_miRNA_expression_boxplot",
  width=12,
  height=10
)

############################################################
# Top20 DE miRNA 单独 ROC/AUC分析
# NAFL vs Control
############################################################


library(pROC)


# 建立真实标签
roc_meta <- meta_f0 %>%
  rownames_to_column("GSM") %>%
  select(
    GSM,
    group
  )


roc_meta$outcome <- ifelse(
  roc_meta$group=="NAFL",
  1,
  0
)



# Top20 miRNA

top20_miRNA <- res_f0_df %>%
  filter(!is.na(padj)) %>%
  filter(!miRNA %in% exclude_miRNAs) %>%
  arrange(padj) %>%
  slice_head(n=20) %>%
  pull(miRNA)



############################################################
# 计算每个miRNA AUC
############################################################


auc_top20 <- map_dfr(
  top20_miRNA,
  function(mir){
    
    expr <- log_expr_f0[
      mir,
      roc_meta$GSM
    ]
    
    
    roc_obj <- roc(
      response = roc_meta$outcome,
      predictor = as.numeric(expr),
      quiet = TRUE
    )
    
    
    tibble(
      miRNA = mir,
      AUC = as.numeric(
        auc(roc_obj)
      ),
      CI_low = as.numeric(
        ci.auc(roc_obj)[1]
      ),
      CI_high = as.numeric(
        ci.auc(roc_obj)[3]
      )
    )
    
  }
) %>%
  arrange(desc(AUC))



print(auc_top20)



write_csv(
  auc_top20,
  file.path(
    tab_dir,
    "Top20_miRNA_individual_AUC.csv"
  )
)
############################################################
# 6. 火山图
############################################################

res_volcano <- res_f0_df %>%
  filter(!miRNA %in% exclude_miRNAs) %>%
  filter(!is.na(padj)) %>%
  mutate(
    change = case_when(
      padj < padj_cut & log2FoldChange >  lfc_cut ~ "Up",
      padj < padj_cut & log2FoldChange < -lfc_cut ~ "Down",
      TRUE ~ "NS"
    ),
    padj_plot = if_else(
      padj == 0,
      min(padj[padj > 0], na.rm = TRUE),
      padj
    ),
    neg_log10_padj = -log10(padj_plot)
  )

n_up   <- sum(res_volcano$change == "Up")
n_down <- sum(res_volcano$change == "Down")
n_ns   <- sum(res_volcano$change == "NS")

top_label <- res_volcano %>%
  filter(change != "NS") %>%
  arrange(padj) %>%
  slice_head(n = 15) %>%
  pull(miRNA)

res_volcano <- res_volcano %>%
  mutate(label = if_else(miRNA %in% top_label, miRNA, NA_character_))

p_volcano <- ggplot(
  res_volcano,
  aes(x = log2FoldChange, y = neg_log10_padj)
) +
  geom_point(aes(color = change), size = 1.4, alpha = 0.65) +
  geom_hline(
    yintercept = -log10(padj_cut),
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.4
  ) +
  geom_vline(
    xintercept = c(-lfc_cut, lfc_cut),
    linetype = "dashed",
    color = "grey40",
    linewidth = 0.4
  ) +
  geom_text_repel(
    aes(label = label),
    size = 3.4,
    fontface = "bold.italic",
    box.padding = 0.5,
    point.padding = 0.3,
    segment.color = "grey40",
    max.overlaps = 50,
    na.rm = TRUE
  ) +
  scale_color_manual(
    values = c(
      "Up" = "#E41A1C",
      "Down" = "#377EB8",
      "NS" = "grey70"
    ),
    labels = c(
      "Up" = paste0("Up (", n_up, ")"),
      "Down" = paste0("Down (", n_down, ")"),
      "NS" = paste0("NS (", n_ns, ")")
    )
  ) +
  labs(
    title = "Volcano Plot: Control vs NAFL",
    subtitle = "miR-3687 excluded",
    x = expression(log[2]~Fold~Change),
    y = expression(-log[10]~adjusted~italic(P)~value),
    color = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = c(0.84, 0.84),
    legend.background = element_rect(fill = alpha("white", 0.8), color = NA),
    panel.grid.minor = element_blank()
  )

print(p_volcano)
save_gg(p_volcano, "Volcano_Control_vs_NAFL_no_miR3687", 8, 6)

############################################################
# 7. 归一化表达矩阵
############################################################

log_expr_f0 <- get_log_norm(counts_f0, meta_f0)

log_expr_f0 <- log_expr_f0[
  !rownames(log_expr_f0) %in% exclude_miRNAs,
  ,
  drop = FALSE
]

write_csv(
  as.data.frame(log_expr_f0) %>%
    rownames_to_column("miRNA"),
  file.path(tab_dir, "log2_normalized_expression_Control_vs_NAFL.csv")
)

############################################################
# 8. 仅针对Control vs NAFL筛选最佳3-miRNA panel
############################################################

panel_dir <- file.path(
  final_dir,
  "Disease_diagnostic_three_miRNA_panel"
)

panel_fig_dir <- file.path(
  panel_dir,
  "figures"
)

panel_tab_dir <- file.path(
  panel_dir,
  "tables"
)

dir.create(
  panel_fig_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  panel_tab_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# 8.1 安全计算AUC
############################################################

safe_auc <- function(response, predictor) {
  
  ok <- complete.cases(
    response,
    predictor
  )
  
  response <- response[ok]
  predictor <- predictor[ok]
  
  if (
    length(response) < 4 ||
    length(unique(response)) < 2 ||
    length(unique(predictor)) < 2
  ) {
    return(NA_real_)
  }
  
  roc_obj <- tryCatch(
    pROC::roc(
      response = response,
      predictor = predictor,
      direction = "auto",
      quiet = TRUE
    ),
    error = function(e) NULL
  )
  
  if (is.null(roc_obj)) {
    return(NA_real_)
  }
  
  as.numeric(
    pROC::auc(roc_obj)
  )
}


############################################################
# 8.2 计算panel内miRNA的最大相关系数
############################################################

max_abs_spearman <- function(
    expr_mat,
    features
) {
  
  features <- intersect(
    features,
    rownames(expr_mat)
  )
  
  if (length(features) < 2) {
    return(0)
  }
  
  cor_mat <- suppressWarnings(
    cor(
      t(
        expr_mat[
          features,
          ,
          drop = FALSE
        ]
      ),
      method = "spearman",
      use = "pairwise.complete.obs"
    )
  )
  
  cor_value <- abs(
    cor_mat[
      upper.tri(cor_mat)
    ]
  )
  
  cor_value <- cor_value[
    is.finite(cor_value)
  ]
  
  if (length(cor_value) == 0) {
    return(NA_real_)
  }
  
  max(cor_value)
}
############################################################
# 计算组内最大绝对Spearman相关系数
# 分别在Control和NAFL内部计算
############################################################

max_abs_spearman_within_group <- function(
    expr_mat,
    meta_df,
    features,
    group_col = "group"
) {
  
  features <- intersect(
    features,
    rownames(expr_mat)
  )
  
  samples <- intersect(
    colnames(expr_mat),
    rownames(meta_df)
  )
  
  if (length(features) < 2) {
    return(0)
  }
  
  meta_use <- meta_df[
    samples,
    ,
    drop = FALSE
  ]
  
  expr_use <- expr_mat[
    features,
    samples,
    drop = FALSE
  ]
  
  group_value <- as.character(
    meta_use[
      ,
      group_col,
      drop = TRUE
    ]
  )
  
  group_levels <- unique(
    group_value[
      !is.na(group_value)
    ]
  )
  
  all_cor_values <- numeric(0)
  
  for (current_group in group_levels) {
    
    group_samples <- rownames(meta_use)[
      group_value == current_group
    ]
    
    if (length(group_samples) < 4) {
      next
    }
    
    cor_matrix <- suppressWarnings(
      cor(
        t(
          expr_use[
            ,
            group_samples,
            drop = FALSE
          ]
        ),
        method = "spearman",
        use = "pairwise.complete.obs"
      )
    )
    
    current_cor_values <- abs(
      cor_matrix[
        upper.tri(cor_matrix)
      ]
    )
    
    current_cor_values <- current_cor_values[
      is.finite(current_cor_values)
    ]
    
    all_cor_values <- c(
      all_cor_values,
      current_cor_values
    )
  }
  
  if (length(all_cor_values) == 0) {
    return(NA_real_)
  }
  
  max(all_cor_values)
}

############################################################
# 8.3 3-miRNA panel重复分层交叉验证函数
############################################################

evaluate_disease_panel_cv <- function(
    features,
    expr_mat = log_expr_f0,
    meta_df = meta_f0,
    k = 5,
    repeats = 30,
    seed = 123
) {
  
  features <- intersect(
    features,
    rownames(expr_mat)
  )
  
  if (length(features) != 3) {
    return(
      tibble(
        mean_CV_AUC = NA_real_,
        sd_CV_AUC = NA_real_,
        median_CV_AUC = NA_real_,
        minimum_CV_AUC = NA_real_,
        maximum_CV_AUC = NA_real_,
        valid_repeats = 0
      )
    )
  }
  
  sample_names <- intersect(
    rownames(meta_df),
    colnames(expr_mat)
  )
  
  model_x <- t(
    expr_mat[
      features,
      sample_names,
      drop = FALSE
    ]
  ) %>%
    as.data.frame()
  
  safe_names <- paste0(
    "x",
    seq_along(features)
  )
  
  colnames(model_x) <- safe_names
  
  model_df <- model_x
  
  model_df$GSM <- rownames(model_x)
  
  model_df$outcome <- as.numeric(
    meta_df[
      model_df$GSM,
      "outcome",
      drop = TRUE
    ]
  )
  
  model_df <- model_df %>%
    select(
      GSM,
      outcome,
      all_of(safe_names)
    )
  
  class_number <- table(
    model_df$outcome
  )
  
  k_use <- min(
    k,
    min(class_number)
  )
  
  if (k_use < 2) {
    return(
      tibble(
        mean_CV_AUC = NA_real_,
        sd_CV_AUC = NA_real_,
        median_CV_AUC = NA_real_,
        minimum_CV_AUC = NA_real_,
        maximum_CV_AUC = NA_real_,
        valid_repeats = 0
      )
    )
  }
  
  model_formula <- as.formula(
    paste(
      "outcome ~",
      paste(
        safe_names,
        collapse = " + "
      )
    )
  )
  
  repeat_auc <- rep(
    NA_real_,
    repeats
  )
  
  for (repeat_i in seq_len(repeats)) {
    
    set.seed(
      seed + repeat_i
    )
    
    folds <- make_stratified_folds(
      model_df$outcome,
      k = k_use
    )
    
    cv_prediction <- rep(
      NA_real_,
      nrow(model_df)
    )
    
    for (fold_i in seq_along(folds)) {
      
      test_index <- folds[[fold_i]]
      
      train_index <- setdiff(
        seq_len(nrow(model_df)),
        test_index
      )
      
      train_df <- model_df[
        train_index,
        ,
        drop = FALSE
      ]
      
      test_df <- model_df[
        test_index,
        ,
        drop = FALSE
      ]
      
      scaled_result <- scale_train_test(
        train_df = train_df,
        test_df = test_df,
        vars = safe_names
      )
      
      train_scaled <- scaled_result$train
      test_scaled  <- scaled_result$test
      
      fit_model <- tryCatch(
        suppressWarnings(
          glm(
            model_formula,
            data = train_scaled,
            family = binomial()
          )
        ),
        error = function(e) NULL
      )
      
      if (is.null(fit_model)) {
        next
      }
      
      test_prediction <- tryCatch(
        suppressWarnings(
          predict(
            fit_model,
            newdata = test_scaled,
            type = "response"
          )
        ),
        error = function(e) {
          rep(
            NA_real_,
            nrow(test_scaled)
          )
        }
      )
      
      cv_prediction[
        test_index
      ] <- test_prediction
    }
    
    repeat_auc[
      repeat_i
    ] <- safe_auc(
      response = model_df$outcome,
      predictor = cv_prediction
    )
  }
  
  valid_auc <- repeat_auc[
    is.finite(repeat_auc)
  ]
  
  if (length(valid_auc) == 0) {
    return(
      tibble(
        mean_CV_AUC = NA_real_,
        sd_CV_AUC = NA_real_,
        median_CV_AUC = NA_real_,
        minimum_CV_AUC = NA_real_,
        maximum_CV_AUC = NA_real_,
        valid_repeats = 0
      )
    )
  }
  
  tibble(
    mean_CV_AUC = mean(valid_auc),
    sd_CV_AUC = sd(valid_auc),
    median_CV_AUC = median(valid_auc),
    minimum_CV_AUC = min(valid_auc),
    maximum_CV_AUC = max(valid_auc),
    valid_repeats = length(valid_auc)
  )
}


############################################################
# 9. 计算全部miRNA的单独AUC
############################################################

disease_outcome <- meta_f0[
  colnames(log_expr_f0),
  "outcome",
  drop = TRUE
]


auc_disease_all <- map_dfr(
  rownames(log_expr_f0),
  function(mir) {
    
    expression_value <- as.numeric(
      log_expr_f0[
        mir,
        colnames(log_expr_f0),
        drop = TRUE
      ]
    )
    
    ok <- complete.cases(
      disease_outcome,
      expression_value
    )
    
    if (
      sum(ok) < 4 ||
      length(unique(disease_outcome[ok])) < 2 ||
      length(unique(expression_value[ok])) < 2
    ) {
      return(
        tibble(
          miRNA = mir,
          AUC = NA_real_,
          CI_low = NA_real_,
          CI_high = NA_real_
        )
      )
    }
    
    roc_obj <- tryCatch(
      pROC::roc(
        response = disease_outcome[ok],
        predictor = expression_value[ok],
        direction = "auto",
        quiet = TRUE
      ),
      error = function(e) NULL
    )
    
    if (is.null(roc_obj)) {
      return(
        tibble(
          miRNA = mir,
          AUC = NA_real_,
          CI_low = NA_real_,
          CI_high = NA_real_
        )
      )
    }
    
    auc_ci <- tryCatch(
      as.numeric(
        pROC::ci.auc(roc_obj)
      ),
      error = function(e) {
        c(
          NA_real_,
          NA_real_,
          NA_real_
        )
      }
    )
    
    tibble(
      miRNA = mir,
      AUC = as.numeric(
        pROC::auc(roc_obj)
      ),
      CI_low = auc_ci[1],
      CI_high = auc_ci[3]
    )
  }
) %>%
  arrange(
    desc(AUC)
  )


write_csv(
  auc_disease_all,
  file.path(
    panel_tab_dir,
    "All_miRNA_individual_AUC_Control_vs_NAFL.csv"
  )
)


cat(
  "\n单miRNA AUC排名前20：\n"
)

print(
  head(
    auc_disease_all,
    20
  )
)

############################################################
# 10. 先筛选“高表达量且单miRNA AUC高”的候选miRNA
#
# 第一步：高原始表达量、较高中位表达量和较高检出率
# 第二步：单miRNA AUC达到high_auc_cut
# 第三步：在通过硬筛的miRNA中，按表达可检出性和AUC综合排序
#
# DESeq2的padj和|log2FC|不再作为首轮硬筛条件，避免优先选到
# 统计显著但原始表达过低、临床上难以稳定检出的miRNA。
# 但仍输出pass_DE_support，供生物学解释和结果复核。
############################################################

if (!exists("expression_summary")) {
  stop(
    paste0(
      "未找到expression_summary。",
      "请从脚本开头完整运行。"
    )
  )
}

required_preselection_columns <- c(
  "miRNA",
  "mean_raw_count",
  "median_raw_count",
  "detection_rate",
  "min_group_detection_rate",
  "expression_score"
)

if (!all(required_preselection_columns %in% colnames(expression_summary))) {
  stop(
    paste0(
      "expression_summary缺少高表达预筛所需列：",
      paste(
        setdiff(
          required_preselection_columns,
          colnames(expression_summary)
        ),
        collapse = ", "
      )
    )
  )
}

# 汇总所有真实miRNA的表达量、AUC和DESeq2证据，作为完整审计表
miRNA_preselection_audit <- res_f0_df %>%
  filter(
    !miRNA %in% exclude_miRNAs
  ) %>%
  inner_join(
    auc_disease_all,
    by = "miRNA"
  ) %>%
  inner_join(
    expression_summary,
    by = "miRNA"
  ) %>%
  mutate(
    pass_mean_raw_count =
      is.finite(mean_raw_count) &
      mean_raw_count >= high_expr_min_mean_raw_count,

    pass_median_raw_count =
      is.finite(median_raw_count) &
      median_raw_count >= high_expr_min_median_raw_count,

    pass_overall_detection =
      is.finite(detection_rate) &
      detection_rate >= high_expr_min_detection_rate,

    pass_group_detection =
      is.finite(min_group_detection_rate) &
      min_group_detection_rate >= high_expr_min_group_detection_rate,

    pass_high_expression =
      pass_mean_raw_count &
      pass_median_raw_count &
      pass_overall_detection &
      pass_group_detection,

    pass_high_AUC =
      is.finite(AUC) &
      AUC >= high_auc_cut,

    pass_DE_support =
      !is.na(padj) &
      padj < panel_padj_cut &
      abs(log2FoldChange) >= panel_lfc_cut,

    pass_expression_and_AUC =
      pass_high_expression & pass_high_AUC,

    expression_direction = case_when(
      log2FoldChange > 0 ~ "Higher in NAFL",
      log2FoldChange < 0 ~ "Lower in NAFL",
      TRUE ~ "No change"
    )
  )

write_csv(
  miRNA_preselection_audit,
  file.path(
    panel_tab_dir,
    "All_miRNA_high_expression_high_AUC_preselection_audit.csv"
  )
)

cat("\n高表达量预筛通过数量：", sum(miRNA_preselection_audit$pass_high_expression), "\n")
cat("高AUC预筛通过数量：", sum(miRNA_preselection_audit$pass_high_AUC), "\n")
cat(
  "同时通过高表达量和高AUC预筛的数量：",
  sum(miRNA_preselection_audit$pass_expression_and_AUC),
  "\n"
)

# 只在通过两项硬筛的miRNA中重新计算排名
# 这样进入组合搜索的候选首先保证可检出性和单标志物诊断能力
high_expression_high_auc_candidates <- miRNA_preselection_audit %>%
  filter(
    pass_expression_and_AUC,
    is.finite(expression_score),
    is.finite(AUC)
  ) %>%
  mutate(
    rank_auc_within_passed = safe_percent_rank(AUC),

    preselection_score =
      preselect_expression_weight * expression_score +
      preselect_auc_weight * rank_auc_within_passed,

    # 为兼容后续旧代码，candidate_score等同于本次预筛综合评分
    candidate_score = preselection_score
  ) %>%
  arrange(
    desc(preselection_score),
    desc(expression_score),
    desc(AUC),
    desc(mean_raw_count),
    desc(detection_rate),
    desc(pass_DE_support),
    padj
  )

# 后续组合筛选继续使用disease_candidates对象名
# 但其含义已改为“高表达量且高AUC的预筛候选”
disease_candidates <- high_expression_high_auc_candidates

write_csv(
  disease_candidates,
  file.path(
    panel_tab_dir,
    "High_expression_high_AUC_miRNA_candidates.csv"
  )
)

# 同时保留旧文件名，避免下游读取脚本中断
write_csv(
  disease_candidates,
  file.path(
    panel_tab_dir,
    "Disease_diagnostic_miRNA_candidates.csv"
  )
)

cat(
  "\n最终进入排序的高表达量且高AUC候选miRNA数量：",
  nrow(disease_candidates),
  "\n"
)

print(
  disease_candidates %>%
    select(
      miRNA,
      mean_raw_count,
      median_raw_count,
      detection_rate,
      detection_rate_Control,
      detection_rate_NAFL,
      min_group_detection_rate,
      expression_cv,
      expression_score,
      AUC,
      CI_low,
      CI_high,
      baseMean,
      log2FoldChange,
      padj,
      pass_DE_support,
      expression_direction,
      preselection_score
    ),
  n = Inf,
  width = Inf
)

if (nrow(disease_candidates) < 3) {
  stop(
    paste0(
      "当前高表达量且高AUC条件下只有",
      nrow(disease_candidates),
      "个候选miRNA，无法构建3-miRNA panel。",
      "建议按以下顺序小幅放宽：",
      "high_expr_min_mean_raw_count从100降至50；",
      "high_expr_min_detection_rate从0.80降至0.70；",
      "最后再把high_auc_cut从0.80降至0.75。"
    )
  )
}


############################################################
# 11. 从“高表达量且高AUC”候选中构建三组合搜索池
############################################################

candidate_number <- min(
  panel_candidate_n,
  nrow(disease_candidates)
)


panel_candidate_pool <- disease_candidates %>%
  slice_head(
    n = candidate_number
  )


candidate_miRNAs <- panel_candidate_pool$miRNA


write_csv(
  panel_candidate_pool,
  file.path(
    panel_tab_dir,
    "High_expression_high_AUC_candidate_pool_for_three_miRNA_panel.csv"
  )
)


cat(
  "\n进入三组合搜索的候选miRNA：\n"
)

print(
  candidate_miRNAs
)


############################################################
# 11.1 输出候选miRNA相关性矩阵
############################################################

candidate_cor_matrix <- cor(
  t(
    log_expr_f0[
      candidate_miRNAs,
      ,
      drop = FALSE
    ]
  ),
  method = "spearman",
  use = "pairwise.complete.obs"
)


write_csv(
  as.data.frame(candidate_cor_matrix) %>%
    rownames_to_column("miRNA"),
  file.path(
    panel_tab_dir,
    "Candidate_miRNA_Spearman_correlation_matrix.csv"
  )
)


############################################################
# 12. 遍历全部3-miRNA组合
#
# 修改说明：
# 1. 不再使用相关性作为硬性淘汰条件
# 2. 所有组合都进入重复交叉验证
# 3. 相关性只记录在结果表中
# 4. 在CV-AUC接近时，优先选择更稳定、相关性更低的组合
############################################################

panel_combinations <- combn(
  candidate_miRNAs,
  3,
  simplify = FALSE
)

cat(
  "\n候选miRNA数量：",
  length(candidate_miRNAs),
  "\n"
)

cat(
  "需要评估的3-miRNA组合数量：",
  length(panel_combinations),
  "\n"
)


############################################################
# 12.1 评估全部三miRNA组合
############################################################

panel_search_results <- map_dfr(
  seq_along(panel_combinations),
  function(combo_i) {
    
    selected_features <- panel_combinations[[combo_i]]
    
    ########################################################
    # 显示计算进度
    ########################################################
    
    if (
      combo_i == 1 ||
      combo_i %% 20 == 0 ||
      combo_i == length(panel_combinations)
    ) {
      
      cat(
        "正在评估组合：",
        combo_i,
        "/",
        length(panel_combinations),
        "\n"
      )
    }
    
    
    ########################################################
    # 计算全部样本混合后的相关性
    ########################################################
    
    pooled_rho <- max_abs_spearman(
      expr_mat = log_expr_f0,
      features = selected_features
    )
    
    
    ########################################################
    # 分别在Control和NAFL内部计算相关性
    ########################################################
    
    within_group_rho <- max_abs_spearman_within_group(
      expr_mat = log_expr_f0,
      meta_df = meta_f0,
      features = selected_features,
      group_col = "group"
    )
    
    
    ########################################################
    # 只记录相关性，不再删除组合
    ########################################################
    
    pass_correlation_085 <- (
      is.finite(within_group_rho) &&
        within_group_rho <= 0.85
    )
    
    pass_correlation_090 <- (
      is.finite(within_group_rho) &&
        within_group_rho <= 0.90
    )
    
    correlation_level <- case_when(
      !is.finite(within_group_rho) ~
        "Correlation unavailable",
      
      within_group_rho >= 0.98 ~
        "Extremely high, rho >= 0.98",
      
      within_group_rho >= 0.90 ~
        "High, rho 0.90-0.98",
      
      within_group_rho >= 0.85 ~
        "Moderately high, rho 0.85-0.90",
      
      TRUE ~
        "Relatively acceptable, rho < 0.85"
    )
    
    
    ########################################################
    # 重复分层交叉验证
    ########################################################
    
    cv_result <- evaluate_disease_panel_cv(
      features = selected_features,
      expr_mat = log_expr_f0,
      meta_df = meta_f0,
      k = k_folds,
      repeats = panel_cv_repeats,
      seed = 10000 + combo_i
    )
    
    
    ########################################################
    # CV无法正常计算时才删除该组合
    ########################################################
    
    if (
      length(cv_result$mean_CV_AUC) == 0 ||
      !is.finite(cv_result$mean_CV_AUC[1])
    ) {
      return(NULL)
    }
    
    
    ########################################################
    # 构建全数据模型
    # 计算表观AUC，仅用于辅助展示
    ########################################################
    
    full_x <- t(
      log_expr_f0[
        selected_features,
        rownames(meta_f0),
        drop = FALSE
      ]
    ) %>%
      as.data.frame()
    
    colnames(full_x) <- c(
      "x1",
      "x2",
      "x3"
    )
    
    full_x$GSM <- rownames(full_x)
    
    full_x$outcome <- meta_f0[
      full_x$GSM,
      "outcome",
      drop = TRUE
    ]
    
    
    ########################################################
    # 全数据标准化
    # 只用于计算表观AUC
    ########################################################
    
    full_x <- full_x %>%
      mutate(
        across(
          c(
            x1,
            x2,
            x3
          ),
          ~ {
            current_sd <- sd(.x, na.rm = TRUE)
            
            if (
              is.na(current_sd) ||
              current_sd == 0
            ) {
              rep(
                0,
                length(.x)
              )
            } else {
              as.numeric(
                scale(.x)
              )
            }
          }
        )
      )
    
    
    ########################################################
    # 全数据Logistic回归
    ########################################################
    
    full_model <- tryCatch(
      suppressWarnings(
        glm(
          outcome ~ x1 + x2 + x3,
          data = full_x,
          family = binomial()
        )
      ),
      error = function(e) NULL
    )
    
    
    ########################################################
    # 表观AUC
    ########################################################
    
    if (is.null(full_model)) {
      
      apparent_auc <- NA_real_
      
    } else {
      
      full_prediction <- tryCatch(
        suppressWarnings(
          predict(
            full_model,
            type = "response"
          )
        ),
        error = function(e) {
          rep(
            NA_real_,
            nrow(full_x)
          )
        }
      )
      
      apparent_auc <- safe_auc(
        response = full_x$outcome,
        predictor = full_prediction
      )
    }
    
    
    ########################################################
    # 汇总当前panel的原始表达丰度与可检出性
    ########################################################

    panel_expression_result <- summarize_panel_expression(
      features = selected_features,
      expression_df = expression_summary
    )

    ########################################################
    # 返回当前组合的结果
    ########################################################
    
    tibble(
      combo_i = combo_i,
      
      miRNA1 = selected_features[1],
      miRNA2 = selected_features[2],
      miRNA3 = selected_features[3],
      
      mean_CV_AUC =
        cv_result$mean_CV_AUC[1],
      
      sd_CV_AUC =
        cv_result$sd_CV_AUC[1],
      
      median_CV_AUC =
        cv_result$median_CV_AUC[1],
      
      minimum_CV_AUC =
        cv_result$minimum_CV_AUC[1],
      
      maximum_CV_AUC =
        cv_result$maximum_CV_AUC[1],
      
      valid_repeats =
        cv_result$valid_repeats[1],
      
      apparent_AUC =
        apparent_auc,

      panel_mean_expression_score =
        panel_expression_result$panel_mean_expression_score[1],

      panel_min_expression_score =
        panel_expression_result$panel_min_expression_score[1],

      panel_mean_raw_count =
        panel_expression_result$panel_mean_raw_count[1],

      panel_median_raw_count =
        panel_expression_result$panel_median_raw_count[1],

      panel_mean_detection_rate =
        panel_expression_result$panel_mean_detection_rate[1],

      panel_min_detection_rate =
        panel_expression_result$panel_min_detection_rate[1],
      
      max_abs_spearman_rho =
        within_group_rho,
      
      pooled_abs_spearman_rho =
        pooled_rho,
      
      pass_correlation_085 =
        pass_correlation_085,
      
      pass_correlation_090 =
        pass_correlation_090,
      
      correlation_level =
        correlation_level
    )
  }
)


############################################################
# 12.2 检查是否存在有效组合
############################################################

cat(
  "\n成功完成交叉验证的组合数量：\n"
)

print(
  nrow(panel_search_results)
)


if (nrow(panel_search_results) == 0) {
  
  stop(
    paste0(
      "所有组合的交叉验证均未得到有效AUC。",
      "此时问题已与相关性无关，",
      "需要检查evaluate_disease_panel_cv函数和样本分组。"
    )
  )
}


############################################################
# 13. 对全部组合排序
#
# 第一标准：平均交叉验证AUC越高越好
# 第二标准：panel平均表达可检出性越高越好
# 第三标准：panel最低检出率越高越好
# 第四标准：交叉验证AUC标准差越小越好
# 第五标准：最低一次交叉验证AUC越高越好
# 第六标准：组内相关性越低越好
############################################################

panel_search_results <- panel_search_results %>%
  arrange(
    desc(mean_CV_AUC),
    desc(panel_mean_expression_score),
    desc(panel_min_detection_rate),
    sd_CV_AUC,
    desc(minimum_CV_AUC),
    max_abs_spearman_rho
  )


write_csv(
  panel_search_results,
  file.path(
    panel_tab_dir,
    "All_three_miRNA_panel_CV_results_no_hard_correlation_filter.csv"
  )
)


cat(
  "\n按照平均CV-AUC排序的前20个组合：\n"
)

print(
  panel_search_results %>%
    slice_head(n = 20),
  n = 20,
  width = Inf
)


############################################################
# 13.1 汇总相关性情况
############################################################

correlation_summary <- panel_search_results %>%
  dplyr::group_by(
    correlation_level
  ) %>%
  dplyr::summarise(
    panel_number = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    correlation_level
  )


write_csv(
  correlation_summary,
  file.path(
    panel_tab_dir,
    "Three_miRNA_panel_correlation_summary.csv"
  )
)


cat(
  "\n全部组合的相关性等级：\n"
)

print(
  correlation_summary
)


############################################################
# 13.2 找到最高平均CV-AUC
############################################################

best_mean_auc <- max(
  panel_search_results$mean_CV_AUC,
  na.rm = TRUE
)


cat(
  "\n所有组合中的最高平均CV-AUC：",
  round(
    best_mean_auc,
    4
  ),
  "\n"
)


############################################################
# 13.3 筛选与最高AUC相差不超过0.01的组合
#
# 这些组合在预测能力上可视为近似。
# 因候选已通过高表达量和高AUC硬筛，在近似最佳组合中进一步优先选择
# panel平均表达可检出性和最低检出率更高者，再考虑CV波动与相关性。
############################################################

near_best_panels <- panel_search_results %>%
  filter(
    mean_CV_AUC >=
      best_mean_auc - 0.01
  ) %>%
  arrange(
    desc(panel_mean_expression_score),
    desc(panel_min_detection_rate),
    sd_CV_AUC,
    max_abs_spearman_rho,
    desc(minimum_CV_AUC),
    desc(mean_CV_AUC)
  )


write_csv(
  near_best_panels,
  file.path(
    panel_tab_dir,
    "Near_best_panels_within_0.01_CV_AUC.csv"
  )
)


cat(
  "\n与最高平均CV-AUC相差不超过0.01的组合数量：",
  nrow(near_best_panels),
  "\n"
)


cat(
  "\n近似最佳组合：\n"
)

print(
  near_best_panels,
  n = Inf,
  width = Inf
)


############################################################
# 14. 提取最终最佳3-miRNA panel
#
# 选择规则：
# 1. 平均CV-AUC处于最高值减0.01以内
# 2. 优先选择panel平均表达可检出性更高的组合
# 3. 再优先选择panel最低检出率更高的组合
# 4. 再优先选择CV标准差更小的组合
# 5. 随后考虑组内相关性和最低一次CV-AUC
############################################################

best_panel <- c(
  near_best_panels$miRNA1[1],
  near_best_panels$miRNA2[1],
  near_best_panels$miRNA3[1]
)


best_panel_summary <- near_best_panels %>%
  slice_head(n = 1) %>%
  mutate(
    panel = paste(
      miRNA1,
      miRNA2,
      miRNA3,
      sep = " + "
    ),
    
    highest_mean_CV_AUC_among_all_panels =
      best_mean_auc,
    
    difference_from_highest_mean_AUC =
      best_mean_auc -
      mean_CV_AUC,
    
    selection_rule = paste0(
      "Mean CV-AUC within 0.01 of the maximum; ",
      "then lowest CV SD, higher raw-expression detectability, ",
      "higher detection rate, lower within-group correlation, ",
      "and higher minimum CV-AUC"
    )
  )


write_csv(
  best_panel_summary,
  file.path(
    panel_tab_dir,
    "BEST_three_miRNA_diagnostic_panel.csv"
  )
)


cat(
  "\n最终选择的3-miRNA疾病诊断panel：\n"
)

print(
  best_panel
)


cat(
  "\n最佳panel详细结果：\n"
)

print(
  best_panel_summary,
  width = Inf
)


############################################################
# 15. 输出最佳panel中单个miRNA结果
############################################################

best_panel_individual_results <- disease_candidates %>%
  filter(
    miRNA %in% best_panel
  ) %>%
  arrange(
    match(
      miRNA,
      best_panel
    )
  )


write_csv(
  best_panel_individual_results,
  file.path(
    panel_tab_dir,
    "BEST_panel_individual_miRNA_results.csv"
  )
)


############################################################
# 16. 最佳panel中3个miRNA的表达箱线图
#     增加Control vs NAFL显著性标注
############################################################

best_panel_box_df <- as.data.frame(
  log_expr_f0[
    best_panel,
    ,
    drop = FALSE
  ]
) %>%
  rownames_to_column("miRNA") %>%
  pivot_longer(
    cols = -miRNA,
    names_to = "GSM",
    values_to = "expression"
  ) %>%
  left_join(
    meta_f0 %>%
      rownames_to_column("GSM") %>%
      dplyr::select(
        GSM,
        group
      ),
    by = "GSM"
  ) %>%
  mutate(
    miRNA = factor(
      miRNA,
      levels = best_panel
    ),
    group = factor(
      group,
      levels = c(
        "Control",
        "NAFL"
      )
    )
  ) %>%
  filter(
    !is.na(group),
    !is.na(expression)
  )


############################################################
# 16.1 对每个miRNA进行Wilcoxon秩和检验
############################################################

panel_significance_df <- best_panel_box_df %>%
  group_by(miRNA) %>%
  group_modify(
    ~ {
      
      test_result <- wilcox.test(
        expression ~ group,
        data = .x,
        exact = FALSE
      )
      
      tibble(
        p = test_result$p.value
      )
    }
  ) %>%
  ungroup() %>%
  mutate(
    # 对3个miRNA的P值进行BH校正
    p.adj = p.adjust(
      p,
      method = "BH"
    ),
    
    # 转换成显著性星号
    p.adj.signif = case_when(
      p.adj < 0.0001 ~ "****",
      p.adj < 0.001  ~ "***",
      p.adj < 0.01   ~ "**",
      p.adj < 0.05   ~ "*",
      TRUE           ~ "ns"
    ),
    
    group1 = "Control",
    group2 = "NAFL"
  )


############################################################
# 16.2 计算每个分面的显著性标注高度
############################################################

panel_y_position <- best_panel_box_df %>%
  group_by(miRNA) %>%
  summarise(
    expression_max = max(
      expression,
      na.rm = TRUE
    ),
    
    expression_min = min(
      expression,
      na.rm = TRUE
    ),
    
    expression_range =
      expression_max -
      expression_min,
    
    y.position =
      expression_max +
      ifelse(
        expression_range > 0,
        expression_range * 0.10,
        0.5
      ),
    
    .groups = "drop"
  )


panel_significance_df <- panel_significance_df %>%
  left_join(
    panel_y_position %>%
      dplyr::select(
        miRNA,
        y.position
      ),
    by = "miRNA"
  )


############################################################
# 16.3 输出显著性检验结果
############################################################

cat(
  "\n最佳3-miRNA panel的组间显著性检验：\n"
)

print(
  panel_significance_df,
  width = Inf
)


write_csv(
  panel_significance_df,
  file.path(
    panel_tab_dir,
    "Best_three_miRNA_expression_significance.csv"
  )
)


############################################################
# 16.4 绘制带显著性标注的箱线图
############################################################

p_best_panel_box <- ggplot(
  best_panel_box_df,
  aes(
    x = group,
    y = expression,
    fill = group
  )
) +
  geom_boxplot(
    width = 0.55,
    outlier.shape = NA,
    linewidth = 0.7,
    alpha = 0.85
  ) +
  geom_jitter(
    aes(
      color = group
    ),
    width = 0.12,
    height = 0,
    size = 1.6,
    alpha = 0.70,
    show.legend = FALSE
  ) +
  
  ##########################################################
# 显著性括号及星号
##########################################################

ggpubr::stat_pvalue_manual(
  panel_significance_df,
  label = "p.adj.signif",
  xmin = "group1",
  xmax = "group2",
  y.position = "y.position",
  tip.length = 0.015,
  bracket.size = 0.6,
  size = 5,
  hide.ns = FALSE
) +
  
  facet_wrap(
    ~ miRNA,
    nrow = 1,
    scales = "free_y"
  ) +
  
  scale_fill_manual(
    values = c(
      "Control" = "#F8766D",
      "NAFL" = "#00BFC4"
    )
  ) +
  
  scale_color_manual(
    values = c(
      "Control" = "#333333",
      "NAFL" = "#333333"
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.03,
        0.17
      )
    )
  ) +
  
  labs(
    title = "Expression of the Best Three-miRNA Diagnostic Panel",
    subtitle = "Control vs NAFL",
    x = NULL,
    y = "log2 normalized expression",
    fill = "Group"
  ) +
  
  theme_bw(
    base_size = 15
  ) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 20
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 16
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 13
    ),
    
    axis.title.y = element_text(
      size = 15
    ),
    
    axis.text = element_text(
      size = 12
    ),
    
    legend.title = element_text(
      size = 14
    ),
    
    legend.text = element_text(
      size = 13
    ),
    
    panel.grid.minor = element_blank()
  )


print(
  p_best_panel_box
)


############################################################
# 16.5 保存图片
############################################################

ggsave(
  filename = file.path(
    panel_fig_dir,
    "Best_three_miRNA_expression_boxplot_with_significance.png"
  ),
  plot = p_best_panel_box,
  width = 10,
  height = 7,
  dpi = 600
)


ggsave(
  filename = file.path(
    panel_fig_dir,
    "Best_three_miRNA_expression_boxplot_with_significance.pdf"
  ),
  plot = p_best_panel_box,
  width = 10,
  height = 7
)
############################################################
# 16.6 最佳panel中3个miRNA的无点小提琴图
############################################################

p_best_panel_violin <- ggplot(
  best_panel_box_df,
  aes(
    x = group,
    y = expression,
    fill = group
  )
) +
  
  ##########################################################
# 小提琴图
##########################################################

geom_violin(
  width = 0.85,
  trim = FALSE,
  scale = "width",
  linewidth = 0.7,
  alpha = 0.75
) +
  
  ##########################################################
# 内嵌箱线图
# 不显示箱线图自身的离群点
##########################################################

geom_boxplot(
  width = 0.16,
  outlier.shape = NA,
  linewidth = 0.65,
  alpha = 0.90
) +
  
  ##########################################################
# 显著性括号及校正后星号
##########################################################

ggpubr::stat_pvalue_manual(
  panel_significance_df,
  label = "p.adj.signif",
  xmin = "group1",
  xmax = "group2",
  y.position = "y.position",
  tip.length = 0.015,
  bracket.size = 0.6,
  size = 5,
  hide.ns = FALSE
) +
  
  facet_wrap(
    ~ miRNA,
    nrow = 1,
    scales = "free_y"
  ) +
  
  scale_fill_manual(
    values = c(
      "Control" = "#F8766D",
      "NAFL" = "#00BFC4"
    )
  ) +
  
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.03,
        0.18
      )
    )
  ) +
  
  labs(
    title = "Expression of the Best Three-miRNA Diagnostic Panel",
    subtitle = "Control vs NAFL",
    x = NULL,
    y = "log2 normalized expression",
    fill = "Group"
  ) +
  
  theme_bw(
    base_size = 15
  ) +
  
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 20
    ),
    
    plot.subtitle = element_text(
      hjust = 0.5,
      size = 16
    ),
    
    strip.text = element_text(
      face = "bold",
      size = 13
    ),
    
    axis.title.y = element_text(
      size = 15
    ),
    
    axis.text = element_text(
      size = 12
    ),
    
    legend.title = element_text(
      size = 14
    ),
    
    legend.text = element_text(
      size = 13
    ),
    
    panel.grid.minor = element_blank()
  )


############################################################
# 显示无点小提琴图
############################################################

print(
  p_best_panel_violin
)


############################################################
# 16.7 保存无点小提琴图
############################################################

ggsave(
  filename = file.path(
    panel_fig_dir,
    "Best_three_miRNA_violin_plot_without_points.png"
  ),
  plot = p_best_panel_violin,
  width = 10,
  height = 7,
  dpi = 600
)


ggsave(
  filename = file.path(
    panel_fig_dir,
    "Best_three_miRNA_violin_plot_without_points.pdf"
  ),
  plot = p_best_panel_violin,
  width = 10,
  height = 7
)
############################################################
# 17. 最佳panel全数据Logistic模型及表观ROC
############################################################

best_model_df <- t(
  log_expr_f0[
    best_panel,
    rownames(meta_f0),
    drop = FALSE
  ]
) %>%
  as.data.frame()


colnames(best_model_df) <- c(
  "x1",
  "x2",
  "x3"
)


best_model_df$GSM <- rownames(best_model_df)
best_model_df$outcome <- meta_f0$outcome


best_model_df <- best_model_df %>%
  mutate(
    across(
      c(
        x1,
        x2,
        x3
      ),
      ~ as.numeric(scale(.x))
    )
  )


best_fit <- glm(
  outcome ~ x1 + x2 + x3,
  data = best_model_df,
  family = binomial()
)


best_model_df$panel_probability <- predict(
  best_fit,
  type = "response"
)


best_roc <- pROC::roc(
  response = best_model_df$outcome,
  predictor = best_model_df$panel_probability,
  quiet = TRUE
)


best_apparent_auc <- as.numeric(
  pROC::auc(best_roc)
)


best_apparent_ci <- as.numeric(
  pROC::ci.auc(best_roc)
)


best_coefficient_table <- summary(
  best_fit
)$coefficients %>%
  as.data.frame() %>%
  rownames_to_column("term") %>%
  mutate(
    miRNA = case_when(
      term == "x1" ~ best_panel[1],
      term == "x2" ~ best_panel[2],
      term == "x3" ~ best_panel[3],
      TRUE ~ term
    )
  )


write_csv(
  best_coefficient_table,
  file.path(
    panel_tab_dir,
    "BEST_panel_logistic_coefficients.csv"
  )
)


write_csv(
  best_model_df,
  file.path(
    panel_tab_dir,
    "BEST_panel_full_data_predictions.csv"
  )
)


best_roc_df <- tibble(
  false_positive_rate =
    1 - best_roc$specificities,
  sensitivity =
    best_roc$sensitivities
)


p_best_roc <- ggplot(
  best_roc_df,
  aes(
    x = false_positive_rate,
    y = sensitivity
  )
) +
  geom_line(
    linewidth = 1.2
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  annotate(
    "text",
    x = 0.62,
    y = 0.18,
    label = paste0(
      "Apparent AUC = ",
      sprintf(
        "%.3f",
        best_apparent_auc
      ),
      "\n95% CI: ",
      sprintf(
        "%.3f",
        best_apparent_ci[1]
      ),
      "-",
      sprintf(
        "%.3f",
        best_apparent_ci[3]
      )
    ),
    hjust = 0,
    size = 4.5
  ) +
  labs(
    title = "ROC of the Best Three-miRNA Diagnostic Panel",
    subtitle = "Control vs NAFL; apparent full-data ROC",
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_bw(
    base_size = 14
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    panel.grid.minor = element_blank()
  )


print(
  p_best_roc
)


ggsave(
  file.path(
    panel_fig_dir,
    "BEST_panel_apparent_ROC.png"
  ),
  p_best_roc,
  width = 7,
  height = 6,
  dpi = 300
)


ggsave(
  file.path(
    panel_fig_dir,
    "BEST_panel_apparent_ROC.pdf"
  ),
  p_best_roc,
  width = 7,
  height = 6
)


############################################################
# 18. 完成提示
############################################################

cat(
  "\n疾病诊断3-miRNA panel筛选完成。\n"
)

cat(
  "最佳panel：",
  paste(
    best_panel,
    collapse = " + "
  ),
  "\n"
)

cat(
  "重复交叉验证平均AUC：",
  round(
    best_panel_summary$mean_CV_AUC,
    3
  ),
  "\n"
)

cat(
  "重复交叉验证AUC标准差：",
  round(
    best_panel_summary$sd_CV_AUC,
    3
  ),
  "\n"
)

cat(
  "全数据表观AUC：",
  round(
    best_apparent_auc,
    3
  ),
  "\n"
)

cat(
  "结果目录：",
  panel_dir,
  "\n"
)