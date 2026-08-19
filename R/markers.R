#' Evaluate individual marker discrimination
#'
#' Computes a fixed-direction rank AUC. For each feature, the reported
#' `discrimination_AUC` is `max(AUC, 1 - AUC)`, while `direction` records
#' whether higher or lower values indicate the positive class in the data used
#' for evaluation.
#'
#' @param expression Numeric feature-by-sample matrix.
#' @param outcome Binary vector aligned to expression columns. Values may be
#'   0/1 or class labels.
#' @param positive Positive-class label when `outcome` is not already 0/1.
#' @return Data frame sorted by decreasing discrimination AUC.
#' @export
evaluate_markers <- function(expression, outcome, positive = NULL) {
  expression <- .as_numeric_matrix(expression, "expression")
  if (is.null(rownames(expression))) {
    stop("Expression must have feature row names.", call. = FALSE)
  }
  if (length(outcome) != ncol(expression)) {
    stop("`outcome` must be aligned to expression columns.", call. = FALSE)
  }
  if (is.numeric(outcome) && all(outcome %in% c(0, 1))) {
    y <- as.integer(outcome)
  } else {
    y <- .check_binary_outcome(outcome, positive)$y
  }
  rows <- lapply(seq_len(nrow(expression)), function(i) {
    auc <- .auc_fixed(y, expression[i, ])
    if (!is.finite(auc)) {
      data.frame(feature = rownames(expression)[i], raw_AUC = NA_real_,
                 discrimination_AUC = NA_real_, direction = NA_character_,
                 stringsAsFactors = FALSE)
    } else {
      data.frame(
        feature = rownames(expression)[i],
        raw_AUC = auc,
        discrimination_AUC = max(auc, 1 - auc),
        direction = if (auc >= 0.5) "higher_in_positive" else "lower_in_positive",
        stringsAsFactors = FALSE
      )
    }
  })
  out <- .bind_rows_base(rows)
  out <- out[order(-out$discrimination_AUC, out$feature), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Preselect detectable and discriminative markers
#'
#' @param detectability Output from [score_detectability()].
#' @param marker_metrics Output from [evaluate_markers()].
#' @param min_mean Minimum mean abundance.
#' @param min_median Minimum median abundance.
#' @param min_detection Minimum overall detection rate.
#' @param min_group_detection Minimum detection rate in the least-detected
#'   group.
#' @param min_auc Minimum discrimination AUC.
#' @param weights Weights for detectability and AUC percentile ranks.
#' @param fallback_n If fewer markers pass all hard filters, optionally return
#'   this many highest-scoring finite markers and flag them as fallback.
#' @return Audit table sorted by decreasing candidate score.
#' @export
preselect_markers <- function(
    detectability,
    marker_metrics,
    min_mean = 100,
    min_median = 20,
    min_detection = 0.80,
    min_group_detection = 0.50,
    min_auc = 0.80,
    weights = c(detectability = 0.50, auc = 0.50),
    fallback_n = 0L
) {
  detectability <- as.data.frame(detectability, stringsAsFactors = FALSE)
  marker_metrics <- as.data.frame(marker_metrics, stringsAsFactors = FALSE)
  needed_d <- c("feature", "mean_value", "median_value", "detection_rate",
                "min_group_detection_rate", "detectability_score")
  needed_m <- c("feature", "discrimination_AUC")
  if (!all(needed_d %in% names(detectability))) {
    stop("`detectability` is missing required columns.", call. = FALSE)
  }
  if (!all(needed_m %in% names(marker_metrics))) {
    stop("`marker_metrics` is missing required columns.", call. = FALSE)
  }
  weights <- .validate_weights(weights, 2L, "weights")
  audit <- merge(detectability, marker_metrics, by = "feature", all = FALSE,
                 sort = FALSE)
  audit$pass_mean <- is.finite(audit$mean_value) & audit$mean_value >= min_mean
  audit$pass_median <- is.finite(audit$median_value) &
    audit$median_value >= min_median
  audit$pass_detection <- is.finite(audit$detection_rate) &
    audit$detection_rate >= min_detection
  audit$pass_group_detection <- is.finite(audit$min_group_detection_rate) &
    audit$min_group_detection_rate >= min_group_detection
  audit$pass_auc <- is.finite(audit$discrimination_AUC) &
    audit$discrimination_AUC >= min_auc
  audit$pass_all <- audit$pass_mean & audit$pass_median &
    audit$pass_detection & audit$pass_group_detection & audit$pass_auc
  audit$auc_rank <- .safe_percent_rank(audit$discrimination_AUC)
  audit$candidate_score <- weights[1L] * audit$detectability_score +
    weights[2L] * audit$auc_rank
  audit$selected_by_fallback <- FALSE

  passed <- audit[audit$pass_all, , drop = FALSE]
  passed <- passed[order(-passed$candidate_score,
                         -passed$detectability_score,
                         -passed$discrimination_AUC), , drop = FALSE]

  fallback_n <- as.integer(fallback_n)[1L]
  if (is.finite(fallback_n) && fallback_n > 0L && nrow(passed) < fallback_n) {
    finite <- audit[is.finite(audit$candidate_score), , drop = FALSE]
    finite <- finite[order(-finite$candidate_score,
                           -finite$detectability_score,
                           -finite$discrimination_AUC), , drop = FALSE]
    extra <- finite[!finite$feature %in% passed$feature, , drop = FALSE]
    extra <- utils::head(extra, max(0L, fallback_n - nrow(passed)))
    if (nrow(extra)) {
      extra$selected_by_fallback <- TRUE
      passed <- rbind(passed, extra)
    }
  }
  rownames(passed) <- NULL
  attr(passed, "audit") <- audit
  passed
}
