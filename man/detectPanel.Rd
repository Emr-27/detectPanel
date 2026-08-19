\name{detectPanel}
\alias{detectPanel}
\alias{discover_panel}
\alias{print.detectPanel_result}
\alias{summary.detectPanel_result}
\alias{plot.detectPanel_result}
\title{Discover, Validate, and Refit a Biomarker Panel}
\usage{
detectPanel(...)

discover_panel(
  counts, metadata, outcome, positive = NULL, sample_id = NULL,
  exclude = character(), exclude_regex = NULL, panel_size = 3L,
  candidate_n = 10L, transform_method = c("log2_cpm", "log2", "none"),
  prior_count = 1, detection_threshold = 10,
  detectability_weights = c(0.5, 0.35, 0.15),
  preselection_weights = c(0.5, 0.5), min_mean = 100,
  min_median = 20, min_detection = 0.8, min_group_detection = 0.5,
  min_auc = 0.8, allow_fallback = FALSE, outer_v = 5L,
  outer_repeats = 2L, inner_v = 5L, inner_repeats = 3L,
  near_best_tolerance = 0.01, seed = 123L, max_combinations = 50000L
)
}
\arguments{
  \item{counts}{Non-negative feature-by-sample matrix.}
  \item{metadata}{Sample metadata data frame.}
  \item{outcome}{Name of the binary outcome column.}
  \item{positive}{Positive-class label.}
  \item{sample_id}{Optional sample-ID column in metadata.}
  \item{exclude}{Exact feature names to remove.}
  \item{exclude_regex}{Optional feature-name regular expression to remove.}
  \item{panel_size}{Number of features per panel.}
  \item{candidate_n}{Maximum candidate-pool size inside each outer split.}
  \item{transform_method}{Modeling transformation.}
  \item{prior_count}{Offset used before log transformation.}
  \item{detection_threshold}{Count threshold defining detection.}
  \item{detectability_weights}{Weights for abundance, detection, and stability.}
  \item{preselection_weights}{Weights for detectability and individual AUC.}
  \item{min_mean, min_median, min_detection, min_group_detection, min_auc}{Training-only candidate filters.}
  \item{allow_fallback}{Whether to fill undersized candidate pools by rank.}
  \item{outer_v, outer_repeats}{Outer validation settings.}
  \item{inner_v, inner_repeats}{Inner panel-search settings.}
  \item{near_best_tolerance}{AUC tolerance before tie breaking.}
  \item{seed}{Random seed.}
  \item{max_combinations}{Maximum combinations per inner search.}
  \item{...}{Arguments passed from \code{detectPanel()} to \code{discover_panel()}.}
}
\value{A \code{detectPanel_result} object containing nested validation and a final fitted model.}
\description{Runs nested validation and refits the most frequently selected exact panel on all samples.}
\examples{
set.seed(1)
y <- rep(c("Control", "Case"), each = 12)
x <- matrix(rpois(10 * 24, 50), 10,
            dimnames = list(paste0("m", 1:10), paste0("s", 1:24)))
x[1:3, y == "Case"] <- x[1:3, y == "Case"] + 50
meta <- data.frame(group = y, row.names = colnames(x))
\dontrun{
fit <- discover_panel(x, meta, "group", positive = "Case",
  candidate_n = 6, min_mean = 5, min_median = 2,
  min_detection = 0.2, min_group_detection = 0.1, min_auc = 0.55,
  outer_v = 3, outer_repeats = 1, inner_v = 3, inner_repeats = 1)
}
}
