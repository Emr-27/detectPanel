#' Export detectPanel results
#'
#' @param x A `detectPanel_result` or `detectPanel_nested` object.
#' @param path Output directory.
#' @param save_model Whether to save the full R object as an RDS file.
#' @return Invisibly returns normalized output path.
#' @export
export_results <- function(x, path, save_model = TRUE) {
  path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  nested <- if (inherits(x, "detectPanel_result")) x$nested else x
  if (!inherits(nested, "detectPanel_nested")) {
    stop("`x` must be a detectPanel result or nested-validation object.",
         call. = FALSE)
  }
  utils::write.csv(nested$outer_summary,
                   file.path(path, "nested_outer_summary.csv"), row.names = FALSE)
  utils::write.csv(nested$split_metrics,
                   file.path(path, "outer_split_metrics.csv"), row.names = FALSE)
  utils::write.csv(nested$predictions,
                   file.path(path, "outer_predictions.csv"), row.names = FALSE)
  utils::write.csv(nested$aggregated_predictions,
                   file.path(path, "sample_level_predictions.csv"),
                   row.names = FALSE)
  utils::write.csv(nested$panel_frequency,
                   file.path(path, "panel_selection_frequency.csv"),
                   row.names = FALSE)
  utils::write.csv(nested$feature_frequency,
                   file.path(path, "feature_selection_frequency.csv"),
                   row.names = FALSE)
  utils::write.csv(nested$split_selections,
                   file.path(path, "outer_split_selections.csv"),
                   row.names = FALSE)
  utils::write.csv(nested$excluded,
                   file.path(path, "excluded_features.csv"), row.names = FALSE)
  if (nrow(nested$errors)) {
    utils::write.csv(nested$errors,
                     file.path(path, "failed_outer_splits.csv"), row.names = FALSE)
  }
  if (inherits(x, "detectPanel_result")) {
    utils::write.csv(x$final_model$coefficient_table,
                     file.path(path, "final_model_coefficients.csv"),
                     row.names = FALSE)
    utils::write.csv(data.frame(feature = x$final_panel),
                     file.path(path, "final_panel.csv"), row.names = FALSE)
  }
  if (isTRUE(save_model)) saveRDS(x, file.path(path, "detectPanel_result.rds"))
  invisible(path)
}
