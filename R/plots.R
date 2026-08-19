#' Plot marker expression by outcome group
#'
#' @param expression Transformed feature-by-sample matrix.
#' @param metadata Sample metadata aligned by row names.
#' @param features Features to display.
#' @param outcome Outcome column used on the x axis.
#' @param sample_id Optional metadata sample-ID column.
#' @return A `ggplot` object.
#' @export
plot_marker_expression <- function(expression, metadata, features, outcome,
                                   sample_id = NULL) {
  expression <- .as_numeric_matrix(expression, "expression")
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (!is.null(sample_id)) {
    if (!sample_id %in% names(metadata)) stop("Unknown `sample_id` column.",
                                              call. = FALSE)
    rownames(metadata) <- as.character(metadata[[sample_id]])
  }
  if (!outcome %in% names(metadata)) stop("Unknown outcome column.",
                                          call. = FALSE)
  samples <- intersect(colnames(expression), rownames(metadata))
  features <- intersect(as.character(features), rownames(expression))
  if (!length(samples) || !length(features)) {
    stop("No shared samples or requested features were found.", call. = FALSE)
  }
  pieces <- lapply(features, function(feature) {
    data.frame(
      sample = samples,
      feature = feature,
      expression = as.numeric(expression[feature, samples]),
      group = metadata[samples, outcome],
      stringsAsFactors = FALSE
    )
  })
  d <- .bind_rows_base(pieces)
  d$feature <- factor(d$feature, levels = features)
  ggplot2::ggplot(d, ggplot2::aes(x = group, y = expression, fill = group)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    ggplot2::geom_jitter(width = 0.15, alpha = 0.7) +
    ggplot2::facet_wrap(~ feature, scales = "free_y") +
    ggplot2::labs(x = NULL, y = "Expression", fill = outcome) +
    ggplot2::theme_bw()
}

#' Plot individual marker ROC curves
#'
#' @param expression Transformed feature-by-sample matrix.
#' @param outcome Binary vector aligned to expression columns.
#' @param features Features to plot.
#' @param positive Positive-class label.
#' @return A `ggplot` object. Each curve is oriented using direction estimated
#'   from the supplied data, so this function is descriptive and should be
#'   called on training or independently locked validation data.
#' @export
plot_marker_roc <- function(expression, outcome, features, positive = NULL) {
  expression <- .as_numeric_matrix(expression, "expression")
  if (length(outcome) != ncol(expression)) stop("Outcome is not aligned.",
                                                call. = FALSE)
  y <- if (is.numeric(outcome) && all(outcome %in% c(0, 1))) {
    as.integer(outcome)
  } else {
    .check_binary_outcome(outcome, positive)$y
  }
  features <- intersect(as.character(features), rownames(expression))
  if (!length(features)) stop("No requested features were found.",
                              call. = FALSE)
  parts <- lapply(features, function(feature) {
    predictor <- as.numeric(expression[feature, ])
    raw_auc <- .auc_fixed(y, predictor)
    oriented <- if (is.finite(raw_auc) && raw_auc < 0.5) -predictor else predictor
    roc <- .roc_points(y, oriented)
    roc$feature <- sprintf("%s (AUC %.3f)", feature,
                           if (is.finite(raw_auc)) max(raw_auc, 1 - raw_auc) else NA_real_)
    roc
  })
  d <- .bind_rows_base(parts)
  ggplot2::ggplot(d, ggplot2::aes(x = fpr, y = sensitivity, color = feature)) +
    ggplot2::geom_line(linewidth = 1) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::coord_equal() +
    ggplot2::labs(x = "1 - Specificity", y = "Sensitivity", color = NULL) +
    ggplot2::theme_bw()
}

#' Plot differential-expression results as a volcano plot
#'
#' @param results Differential-expression result data frame.
#' @param feature_col,lfc_col,padj_col Column names.
#' @param padj_cut Adjusted P-value cutoff.
#' @param lfc_cut Absolute log2-fold-change cutoff.
#' @return A `ggplot` object.
#' @export
plot_volcano <- function(results, feature_col = "feature",
                         lfc_col = "log2FoldChange", padj_col = "padj",
                         padj_cut = 0.05, lfc_cut = 1) {
  results <- as.data.frame(results, stringsAsFactors = FALSE)
  needed <- c(feature_col, lfc_col, padj_col)
  if (!all(needed %in% names(results))) {
    stop("Differential-expression results are missing required columns.",
         call. = FALSE)
  }
  d <- data.frame(
    feature = results[[feature_col]],
    log2FoldChange = as.numeric(results[[lfc_col]]),
    padj = as.numeric(results[[padj_col]]),
    stringsAsFactors = FALSE
  )
  positive <- d$padj[is.finite(d$padj) & d$padj > 0]
  floor_value <- if (length(positive)) min(positive) else .Machine$double.xmin
  d$padj_plot <- ifelse(is.finite(d$padj),
                         pmax(d$padj, floor_value), NA_real_)
  d$neg_log10_padj <- -log10(d$padj_plot)
  d$change <- "NS"
  d$change[is.finite(d$padj) & d$padj < padj_cut &
             d$log2FoldChange > lfc_cut] <- "Up"
  d$change[is.finite(d$padj) & d$padj < padj_cut &
             d$log2FoldChange < -lfc_cut] <- "Down"
  ggplot2::ggplot(d, ggplot2::aes(x = log2FoldChange,
                                  y = neg_log10_padj, color = change)) +
    ggplot2::geom_point(alpha = 0.7) +
    ggplot2::geom_hline(yintercept = -log10(padj_cut), linetype = "dashed") +
    ggplot2::geom_vline(xintercept = c(-lfc_cut, lfc_cut),
                        linetype = "dashed") +
    ggplot2::labs(x = "log2 fold change", y = "-log10 adjusted P value",
                  color = NULL) +
    ggplot2::theme_bw()
}

#' Plot a panel heatmap with optional pheatmap support
#'
#' @param expression Transformed feature-by-sample matrix.
#' @param features Features shown as rows.
#' @param metadata Optional metadata aligned by row names.
#' @param annotation_cols Optional metadata columns used as sample annotations.
#' @param scale Passed to `pheatmap::pheatmap()`.
#' @param ... Additional arguments passed to `pheatmap::pheatmap()`.
#' @return The object returned by `pheatmap::pheatmap()`.
#' @export
plot_panel_heatmap <- function(expression, features, metadata = NULL,
                               annotation_cols = NULL,
                               scale = "row", ...) {
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    stop("Package 'pheatmap' is required for `plot_panel_heatmap()`.",
         call. = FALSE)
  }
  expression <- .as_numeric_matrix(expression, "expression")
  features <- intersect(as.character(features), rownames(expression))
  if (!length(features)) stop("No requested features were found.",
                              call. = FALSE)
  mat <- expression[features, , drop = FALSE]
  annotation <- NULL
  if (!is.null(metadata)) {
    metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
    samples <- intersect(colnames(mat), rownames(metadata))
    if (!length(samples)) stop("Expression and metadata have no shared samples.",
                               call. = FALSE)
    mat <- mat[, samples, drop = FALSE]
    if (!is.null(annotation_cols)) {
      if (!all(annotation_cols %in% names(metadata))) {
        stop("Some `annotation_cols` were not found in metadata.",
             call. = FALSE)
      }
      annotation <- metadata[samples, annotation_cols, drop = FALSE]
    }
  }
  pheatmap::pheatmap(mat, annotation_col = annotation, scale = scale, ...)
}

#' Compute confusion-matrix metrics
#'
#' @param response Binary 0/1 outcome.
#' @param probability Predicted probability for class 1.
#' @param threshold Probability classification threshold.
#' @return One-row data frame with AUC, sensitivity, specificity, accuracy, and
#'   confusion counts.
#' @export
confusion_metrics <- function(response, probability, threshold = 0.5) {
  if (!all(response %in% c(0, 1)) || anyNA(response)) {
    stop("`response` must be a non-missing 0/1 vector.", call. = FALSE)
  }
  .classification_metrics(response, probability, threshold)
}
