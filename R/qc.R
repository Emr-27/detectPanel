#' Flag potential PCA outliers without deleting samples
#'
#' This function is intentionally non-destructive. It returns PCA coordinates
#' and a robust multivariate distance so that suspected samples can be reviewed
#' alongside independent laboratory and sequencing quality-control evidence.
#'
#' @param expression Transformed feature-by-sample matrix.
#' @param group Optional vector aligned to samples, returned for plotting.
#' @param n_components Maximum number of principal components used in the
#'   robust distance.
#' @param threshold Optional robust-distance threshold. Set to `NULL` to return
#'   distances without flagging any sample.
#' @return Data frame with sample names, PCA coordinates, robust distance, and
#'   a logical flag.
#' @export
flag_pca_outliers <- function(expression, group = NULL, n_components = 5L,
                              threshold = NULL) {
  expression <- .as_numeric_matrix(expression, "expression")
  if (is.null(colnames(expression))) {
    stop("Expression must have sample column names.", call. = FALSE)
  }
  if (!is.null(group) && length(group) != ncol(expression)) {
    stop("`group` must be aligned to expression columns.", call. = FALSE)
  }
  finite_rows <- apply(expression, 1L, function(x) all(is.finite(x)))
  variable_rows <- apply(expression, 1L, stats::sd, na.rm = TRUE) > 0
  keep <- finite_rows & variable_rows
  if (sum(keep) < 2L) {
    stop("At least two finite, variable features are required for PCA.",
         call. = FALSE)
  }
  pca <- stats::prcomp(t(expression[keep, , drop = FALSE]),
                       center = TRUE, scale. = TRUE)
  k <- min(as.integer(n_components)[1L], ncol(pca$x), nrow(pca$x) - 1L)
  if (!is.finite(k) || k < 1L) stop("No usable principal components.",
                                    call. = FALSE)
  scores <- pca$x[, seq_len(k), drop = FALSE]
  med <- apply(scores, 2L, stats::median, na.rm = TRUE)
  spread <- apply(scores, 2L, stats::mad, constant = 1.4826, na.rm = TRUE)
  spread[!is.finite(spread) | spread == 0] <- 1
  robust_z <- sweep(sweep(scores, 2L, med, "-"), 2L, spread, "/")
  distance <- sqrt(rowSums(robust_z^2))
  flagged <- if (is.null(threshold)) {
    rep(FALSE, length(distance))
  } else {
    threshold <- as.numeric(threshold)[1L]
    if (!is.finite(threshold) || threshold <= 0) {
      stop("`threshold` must be NULL or a positive finite number.",
           call. = FALSE)
    }
    distance > threshold
  }
  out <- data.frame(
    sample = rownames(scores),
    robust_pca_distance = distance,
    flagged = flagged,
    stringsAsFactors = FALSE
  )
  if (!is.null(group)) out$group <- group
  for (j in seq_len(k)) out[[paste0("PC", j)]] <- scores[, j]
  attr(out, "percent_variance") <- 100 * pca$sdev[seq_len(k)]^2 /
    sum(pca$sdev^2)
  out[order(-out$robust_pca_distance), , drop = FALSE]
}
