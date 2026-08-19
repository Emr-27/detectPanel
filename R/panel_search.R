.max_abs_spearman <- function(expression, features, outcome = NULL,
                              within_group = TRUE) {
  features <- intersect(features, rownames(expression))
  if (length(features) < 2L) return(0)
  groups <- if (is.null(outcome) || !within_group) {
    list(all = seq_len(ncol(expression)))
  } else {
    split(seq_along(outcome), outcome)
  }
  values <- numeric()
  for (idx in groups) {
    if (length(idx) < 4L) next
    cm <- suppressWarnings(stats::cor(
      t(expression[features, idx, drop = FALSE]),
      method = "spearman", use = "pairwise.complete.obs"
    ))
    v <- abs(cm[upper.tri(cm)])
    values <- c(values, v[is.finite(v)])
  }
  if (!length(values)) NA_real_ else max(values)
}

.panel_candidate_summary <- function(features, candidate_table) {
  defaults <- data.frame(
    panel_mean_candidate_score = NA_real_,
    panel_mean_detectability_score = NA_real_,
    panel_min_detection_rate = NA_real_,
    panel_mean_value = NA_real_
  )
  if (is.null(candidate_table)) return(defaults)
  candidate_table <- as.data.frame(candidate_table, stringsAsFactors = FALSE)
  if (!"feature" %in% names(candidate_table)) return(defaults)
  x <- candidate_table[candidate_table$feature %in% features, , drop = FALSE]
  if (nrow(x) != length(features)) return(defaults)
  value <- function(name, fun = mean) {
    if (!name %in% names(x)) return(NA_real_)
    z <- x[[name]]
    z <- z[is.finite(z)]
    if (!length(z)) NA_real_ else fun(z)
  }
  data.frame(
    panel_mean_candidate_score = value("candidate_score"),
    panel_mean_detectability_score = value("detectability_score"),
    panel_min_detection_rate = value("detection_rate", min),
    panel_mean_value = value("mean_value")
  )
}

#' Search small biomarker panels with shared cross-validation splits
#'
#' @param expression Transformed feature-by-sample expression matrix.
#' @param outcome Binary 0/1 vector aligned to expression columns.
#' @param candidates Candidate feature names.
#' @param panel_size Number of features per panel.
#' @param splits Repeated stratified splits from
#'   [make_repeated_stratified_splits()]. If omitted, they are generated once
#'   and reused for every panel.
#' @param candidate_table Optional marker table returned by
#'   [preselect_markers()] for expression-aware tie breaking.
#' @param v Number of folds when `splits` is omitted.
#' @param repeats Number of repeats when `splits` is omitted.
#' @param seed Random seed when `splits` is omitted.
#' @param near_best_tolerance Panels within this amount of the maximum mean
#'   CV AUC are treated as near-best before tie breaking.
#' @param max_combinations Safety limit for exhaustive combinations.
#' @return A `detectPanel_search` object.
#' @export
search_panels <- function(
    expression,
    outcome,
    candidates,
    panel_size = 3L,
    splits = NULL,
    candidate_table = NULL,
    v = 5L,
    repeats = 10L,
    seed = 123L,
    near_best_tolerance = 0.01,
    max_combinations = 50000L
) {
  expression <- .as_numeric_matrix(expression, "expression")
  if (is.null(rownames(expression))) {
    stop("Expression must have feature row names.", call. = FALSE)
  }
  if (length(outcome) != ncol(expression) ||
      !all(outcome %in% c(0, 1)) || anyNA(outcome)) {
    stop("`outcome` must be a non-missing 0/1 vector aligned to columns.",
         call. = FALSE)
  }
  candidates <- unique(as.character(candidates))
  candidates <- candidates[candidates %in% rownames(expression)]
  panel_size <- as.integer(panel_size)[1L]
  if (panel_size < 1L || length(candidates) < panel_size) {
    stop("Not enough candidate features for the requested panel size.",
         call. = FALSE)
  }
  n_combo <- choose(length(candidates), panel_size)
  if (!is.finite(n_combo) || n_combo > max_combinations) {
    stop(sprintf(
      "Panel search would evaluate %.0f combinations; reduce candidates or increase `max_combinations`.",
      n_combo
    ), call. = FALSE)
  }
  if (is.null(splits)) {
    splits <- make_repeated_stratified_splits(outcome, v = v,
                                               repeats = repeats, seed = seed)
  }
  combinations <- utils::combn(candidates, panel_size, simplify = FALSE)
  rows <- vector("list", length(combinations))

  for (i in seq_along(combinations)) {
    features <- combinations[[i]]
    cv <- .evaluate_panel_on_splits(expression, outcome, features, splits)
    candidate_summary <- .panel_candidate_summary(features, candidate_table)
    rows[[i]] <- cbind(
      data.frame(
        combination_id = i,
        panel = .panel_key(features),
        stringsAsFactors = FALSE
      ),
      as.data.frame(as.list(stats::setNames(features,
        paste0("feature", seq_along(features)))), stringsAsFactors = FALSE),
      cv,
      candidate_summary,
      data.frame(
        max_abs_spearman_within_group = .max_abs_spearman(
          expression, features, outcome, within_group = TRUE
        ),
        max_abs_spearman_pooled = .max_abs_spearman(
          expression, features, outcome, within_group = FALSE
        )
      )
    )
  }
  results <- .bind_rows_base(rows)
  results <- results[is.finite(results$mean_CV_AUC), , drop = FALSE]
  if (!nrow(results)) {
    stop("No panel produced a valid cross-validated AUC.", call. = FALSE)
  }

  # Primary ranking. Candidate summaries are used only after predictive
  # performance is declared near-equivalent.
  best_auc <- max(results$mean_CV_AUC, na.rm = TRUE)
  near <- results$mean_CV_AUC >= best_auc - near_best_tolerance
  results$near_best <- near
  results$auc_difference_from_best <- best_auc - results$mean_CV_AUC

  near_results <- results[near, , drop = FALSE]
  order_near <- order(
    -replace(near_results$panel_mean_candidate_score,
             !is.finite(near_results$panel_mean_candidate_score), -Inf),
    -replace(near_results$panel_min_detection_rate,
             !is.finite(near_results$panel_min_detection_rate), -Inf),
    near_results$sd_CV_AUC,
    -near_results$minimum_CV_AUC,
    replace(near_results$max_abs_spearman_within_group,
            !is.finite(near_results$max_abs_spearman_within_group), Inf),
    -near_results$mean_CV_AUC
  )
  near_results <- near_results[order_near, , drop = FALSE]
  selected_panel <- strsplit(near_results$panel[1L], " \\+ ")[[1L]]

  global_order <- order(-results$mean_CV_AUC, results$sd_CV_AUC,
                        -results$minimum_CV_AUC)
  results <- results[global_order, , drop = FALSE]
  rownames(results) <- NULL
  rownames(near_results) <- NULL

  out <- list(
    results = results,
    near_best = near_results,
    best_panel = selected_panel,
    best_auc = best_auc,
    selection_rule = paste(
      "Mean CV AUC within tolerance of maximum; then higher mean candidate",
      "score, higher minimum detection rate, lower CV AUC SD, higher minimum",
      "CV AUC, and lower within-group absolute Spearman correlation."
    ),
    splits = splits,
    settings = list(
      panel_size = panel_size,
      candidate_count = length(candidates),
      combination_count = length(combinations),
      near_best_tolerance = near_best_tolerance
    )
  )
  class(out) <- "detectPanel_search"
  out
}

#' @export
print.detectPanel_search <- function(x, ...) {
  cat("detectPanel panel search\n")
  cat("  Candidates:", x$settings$candidate_count, "\n")
  cat("  Combinations:", x$settings$combination_count, "\n")
  cat("  Best panel:", paste(x$best_panel, collapse = " + "), "\n")
  cat("  Maximum mean CV AUC:", sprintf("%.3f", x$best_auc), "\n")
  invisible(x)
}
