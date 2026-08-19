#' Discover, validate, and refit a biomarker panel
#'
#' Runs [nested_validate_panels()] and then refits the most frequently selected
#' exact panel using all available samples. The all-data model is intended for
#' deployment after the nested validation performance has been reported; its
#' training performance is not an independent estimate of accuracy.
#'
#' @inheritParams nested_validate_panels
#' @return A `detectPanel_result` object containing nested validation and a
#'   final fitted model.
#' @export
discover_panel <- function(
    counts,
    metadata,
    outcome,
    positive = NULL,
    sample_id = NULL,
    exclude = character(),
    exclude_regex = NULL,
    panel_size = 3L,
    candidate_n = 10L,
    transform_method = c("log2_cpm", "log2", "none"),
    prior_count = 1,
    detection_threshold = 10,
    detectability_weights = c(0.50, 0.35, 0.15),
    preselection_weights = c(0.50, 0.50),
    min_mean = 100,
    min_median = 20,
    min_detection = 0.80,
    min_group_detection = 0.50,
    min_auc = 0.80,
    allow_fallback = FALSE,
    outer_v = 5L,
    outer_repeats = 2L,
    inner_v = 5L,
    inner_repeats = 3L,
    near_best_tolerance = 0.01,
    seed = 123L,
    max_combinations = 50000L
) {
  transform_method <- match.arg(transform_method)
  nested <- nested_validate_panels(
    counts = counts,
    metadata = metadata,
    outcome = outcome,
    positive = positive,
    sample_id = sample_id,
    exclude = exclude,
    exclude_regex = exclude_regex,
    panel_size = panel_size,
    candidate_n = candidate_n,
    transform_method = transform_method,
    prior_count = prior_count,
    detection_threshold = detection_threshold,
    detectability_weights = detectability_weights,
    preselection_weights = preselection_weights,
    min_mean = min_mean,
    min_median = min_median,
    min_detection = min_detection,
    min_group_detection = min_group_detection,
    min_auc = min_auc,
    allow_fallback = allow_fallback,
    outer_v = outer_v,
    outer_repeats = outer_repeats,
    inner_v = inner_v,
    inner_repeats = inner_repeats,
    near_best_tolerance = near_best_tolerance,
    seed = seed,
    max_combinations = max_combinations
  )
  validated <- validate_assay_data(
    counts, metadata, outcome, positive = positive, sample_id = sample_id
  )
  filtered <- exclude_features(
    validated$assay, features = exclude, regex = exclude_regex
  )
  final_panel <- strsplit(nested$panel_frequency$panel[1L], " \\+ ")[[1L]]
  final_model <- fit_panel(
    assay = filtered$assay,
    outcome = validated$metadata[[validated$outcome_name]],
    features = final_panel,
    positive = validated$positive,
    transform_method = transform_method,
    prior_count = prior_count,
    threshold_method = "youden"
  )
  out <- list(
    nested = nested,
    final_panel = final_panel,
    final_model = final_model,
    note = paste(
      "Nested out-of-fold metrics estimate internal validation performance.",
      "The final model was refit on all samples and requires external validation."
    )
  )
  class(out) <- "detectPanel_result"
  out
}

#' Main detectPanel workflow
#'
#' This is an alias for [discover_panel()].
#'
#' @inheritParams discover_panel
#' @return A `detectPanel_result` object.
#' @export
detectPanel <- function(...) discover_panel(...)

#' @export
print.detectPanel_result <- function(x, ...) {
  cat("detectPanel discovery result\n")
  cat("  Final panel:", paste(x$final_panel, collapse = " + "), "\n")
  cat("  Mean nested outer AUC:",
      sprintf("%.3f", x$nested$outer_summary$mean_outer_AUC), "\n")
  cat("  Pooled nested AUC:",
      sprintf("%.3f", x$nested$outer_summary$pooled_sample_level_AUC), "\n")
  cat("  Note:", x$note, "\n")
  invisible(x)
}

#' @export
summary.detectPanel_result <- function(object, ...) {
  list(
    final_panel = object$final_panel,
    nested_validation = summary(object$nested),
    final_model = summary(object$final_model),
    note = object$note
  )
}

#' @export
plot.detectPanel_result <- function(x, ...) plot(x$nested, ...)
