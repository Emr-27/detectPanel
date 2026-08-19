\name{score_detectability}
\alias{score_detectability}
\alias{evaluate_markers}
\alias{preselect_markers}
\title{Score and Preselect Candidate Biomarkers}
\usage{
score_detectability(counts, group = NULL, detection_threshold = 10,
  weights = c(abundance = 0.5, detection = 0.35, stability = 0.15),
  stability_transform = c("log2", "none"))

evaluate_markers(expression, outcome, positive = NULL)

preselect_markers(detectability, marker_metrics, min_mean = 100,
  min_median = 20, min_detection = 0.8, min_group_detection = 0.5,
  min_auc = 0.8, weights = c(detectability = 0.5, auc = 0.5),
  fallback_n = 0L)
}
\arguments{
  \item{counts}{Non-negative feature-by-sample count matrix.}
  \item{group}{Optional grouping vector.}
  \item{detection_threshold}{Value defining detection.}
  \item{weights}{Non-negative weights summing to one.}
  \item{stability_transform}{Transformation before robust dispersion scoring.}
  \item{expression}{Transformed feature-by-sample matrix.}
  \item{outcome}{Binary outcome aligned to columns.}
  \item{positive}{Positive-class label.}
  \item{detectability}{Output from \code{score_detectability()}.}
  \item{marker_metrics}{Output from \code{evaluate_markers()}.}
  \item{min_mean, min_median, min_detection, min_group_detection, min_auc}{Candidate filters.}
  \item{fallback_n}{Optional number of ranked markers returned when hard filters are insufficient.}
}
\value{A feature-level data frame.}
\description{Combines abundance, detection, robust stability, and fixed-direction rank AUC.}
