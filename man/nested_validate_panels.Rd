\name{nested_validate_panels}
\alias{nested_validate_panels}
\alias{print.detectPanel_nested}
\alias{summary.detectPanel_nested}
\alias{plot.detectPanel_nested}
\title{Nested Validation of Biomarker Panel Discovery}
\usage{
nested_validate_panels(counts, metadata, outcome, positive = NULL,
  sample_id = NULL, exclude = character(), exclude_regex = NULL,
  panel_size = 3L, candidate_n = 10L,
  transform_method = c("log2_cpm", "log2", "none"), prior_count = 1,
  detection_threshold = 10, detectability_weights = c(0.5, 0.35, 0.15),
  preselection_weights = c(0.5, 0.5), min_mean = 100,
  min_median = 20, min_detection = 0.8, min_group_detection = 0.5,
  min_auc = 0.8, allow_fallback = FALSE, outer_v = 5L,
  outer_repeats = 2L, inner_v = 5L, inner_repeats = 3L,
  near_best_tolerance = 0.01, seed = 123L, max_combinations = 50000L)
}
\arguments{
  \item{counts, metadata, outcome, positive, sample_id, exclude, exclude_regex}{See \code{discover_panel()}.}
  \item{panel_size, candidate_n, transform_method, prior_count}{Model and search settings.}
  \item{detection_threshold, detectability_weights, preselection_weights}{Marker-scoring settings.}
  \item{min_mean, min_median, min_detection, min_group_detection, min_auc}{Training-only candidate filters.}
  \item{allow_fallback}{Whether to fill undersized candidate pools by rank.}
  \item{outer_v, outer_repeats, inner_v, inner_repeats}{Nested resampling settings.}
  \item{near_best_tolerance, seed, max_combinations}{Search controls.}
}
\value{A \code{detectPanel_nested} object with predictions, metrics, selection frequencies, audits, and errors.}
\description{Repeats candidate selection and panel search within every outer training split.
If every outer split fails, the error message includes each split identifier and
its original stage-specific failure message.}
