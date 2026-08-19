#' Score biomarker detectability and within-group stability
#'
#' @param counts Non-negative feature-by-sample count matrix.
#' @param group Optional binary or multi-level grouping vector aligned to assay
#'   columns. Group information is used only for within-group stability and
#'   minimum group detection rate.
#' @param detection_threshold Count threshold defining detection.
#' @param weights Three weights for abundance, detection, and stability.
#' @param stability_transform Transformation used before robust within-group
#'   median absolute deviation: `"log2"` or `"none"`.
#' @return Data frame sorted by decreasing `detectability_score`.
#' @export
score_detectability <- function(
    counts,
    group = NULL,
    detection_threshold = 10,
    weights = c(abundance = 0.50, detection = 0.35, stability = 0.15),
    stability_transform = c("log2", "none")
) {
  counts <- .as_numeric_matrix(counts, "counts")
  if (is.null(rownames(counts))) {
    stop("Counts must have feature row names.", call. = FALSE)
  }
  if (any(counts < 0, na.rm = TRUE)) {
    stop("Counts must be non-negative.", call. = FALSE)
  }
  if (!is.numeric(detection_threshold) || length(detection_threshold) != 1L ||
      !is.finite(detection_threshold) || detection_threshold < 0) {
    stop("`detection_threshold` must be a finite non-negative number.",
         call. = FALSE)
  }
  weights <- .validate_weights(weights, 3L, "weights")
  stability_transform <- match.arg(stability_transform)

  mean_value <- rowMeans(counts, na.rm = TRUE)
  median_value <- apply(counts, 1L, stats::median, na.rm = TRUE)
  detection_rate <- rowMeans(counts >= detection_threshold, na.rm = TRUE)

  if (is.null(group)) {
    group <- factor(rep("all", ncol(counts)))
  } else {
    if (length(group) != ncol(counts)) {
      stop("`group` must be aligned to assay columns.", call. = FALSE)
    }
    group <- factor(group)
    if (anyNA(group)) stop("`group` contains missing values.", call. = FALSE)
  }

  stability_input <- if (stability_transform == "log2") log2(counts + 1) else counts
  group_levels <- levels(group)
  group_detection <- matrix(NA_real_, nrow(counts), length(group_levels),
                            dimnames = list(rownames(counts), group_levels))
  group_dispersion <- group_detection

  for (j in seq_along(group_levels)) {
    idx <- which(group == group_levels[j])
    group_detection[, j] <- rowMeans(
      counts[, idx, drop = FALSE] >= detection_threshold,
      na.rm = TRUE
    )
    group_dispersion[, j] <- apply(
      stability_input[, idx, drop = FALSE], 1L, .safe_mad
    )
  }

  mean_within_group_mad <- apply(group_dispersion, 1L, .safe_mean)
  min_group_detection <- apply(group_detection, 1L, .safe_min)
  abundance_rank <- .safe_percent_rank(log10(mean_value + 1))
  detection_rank <- .safe_percent_rank(detection_rate)
  stability_rank <- .safe_percent_rank(-mean_within_group_mad)
  score <- weights[1L] * abundance_rank +
    weights[2L] * detection_rank + weights[3L] * stability_rank

  out <- data.frame(
    feature = rownames(counts),
    mean_value = mean_value,
    median_value = median_value,
    detection_rate = detection_rate,
    min_group_detection_rate = min_group_detection,
    mean_within_group_mad = mean_within_group_mad,
    abundance_rank = abundance_rank,
    detection_rank = detection_rank,
    stability_rank = stability_rank,
    detectability_score = score,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  for (j in seq_along(group_levels)) {
    safe_name <- make.names(group_levels[j])
    out[[paste0("detection_rate_", safe_name)]] <- group_detection[, j]
    out[[paste0("mad_", safe_name)]] <- group_dispersion[, j]
  }
  out <- out[order(-out$detectability_score, -out$mean_value,
                   -out$detection_rate), , drop = FALSE]
  rownames(out) <- NULL
  out
}
