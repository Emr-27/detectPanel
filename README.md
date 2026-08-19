# detectPanel

`detectPanel` is an R package for leakage-aware discovery of small biomarker
panels from feature-by-sample count matrices.

The package was refactored from a GSE185062 miRNA workflow. Dataset-specific
rules such as `Control`/`NAFL`, removal of `miR-3687`, Windows paths, and fixed
output filenames are no longer embedded in core functions.

## What changed from the original workflow

- Candidate filtering and panel selection are repeated inside each outer
  training split.
- Every candidate panel is compared on the same inner resampling splits.
- `flag_pca_outliers()` provides non-destructive PCA QC flags; samples are not
  automatically removed by a percentile rule.
- AUC direction for panel probabilities is fixed: larger probabilities mean
  the positive class.
- Classification thresholds are selected on training data and applied to outer
  assessment data.
- Fitted models retain transformation, centering, scaling, coefficients, and
  threshold information for prediction on new samples.
- Ordinary logistic regression is used when stable; non-convergent or
  separation-prone fits automatically fall back to an internal L2-penalized
  logistic fit, and the fitted object records which method was used.
- Raw-count detectability scoring is separated from transformed expression used
  for modeling.

## Installation from the source directory

```r
install.packages("detectPanel_0.1.0.tar.gz", repos = NULL, type = "source")
```

During development:

```r
install.packages("remotes")
remotes::install_local("detectPanel")
```

## Minimal example

```r
library(detectPanel)

set.seed(42)
n <- 48
p <- 18
y <- rep(c("Control", "Case"), each = n / 2)

counts <- matrix(
  rpois(p * n, lambda = 80),
  nrow = p,
  dimnames = list(paste0("marker", seq_len(p)), paste0("S", seq_len(n)))
)
counts[1:3, y == "Case"] <- counts[1:3, y == "Case"] + 80
metadata <- data.frame(group = y, row.names = colnames(counts))

result <- discover_panel(
  counts = counts,
  metadata = metadata,
  outcome = "group",
  positive = "Case",
  panel_size = 3,
  candidate_n = 8,
  detection_threshold = 5,
  min_mean = 10,
  min_median = 5,
  min_detection = 0.50,
  min_group_detection = 0.30,
  min_auc = 0.55,
  outer_v = 3,
  outer_repeats = 1,
  inner_v = 3,
  inner_repeats = 1,
  seed = 42
)

result
plot(result, type = "roc")
plot(result, type = "feature_frequency")

# New count matrix: same features in rows, new samples in columns
new_probability <- predict(result, counts[, 1:4], type = "response")
new_class <- predict(result, counts[, 1:4], type = "class")
```

## GSE185062-style call

```r
result <- discover_panel(
  counts = counts_raw,
  metadata = meta,
  outcome = "group",
  positive = "NAFL",
  exclude = "miR-3687",
  exclude_regex = "^(CTRL_|HK_)",
  panel_size = 3,
  candidate_n = 12,
  detection_threshold = 10,
  min_mean = 100,
  min_median = 20,
  min_detection = 0.80,
  min_group_detection = 0.50,
  min_auc = 0.80,
  outer_v = 5,
  outer_repeats = 5,
  inner_v = 5,
  inner_repeats = 10,
  seed = 123
)

export_results(result, "detectPanel-results")
```

Strict thresholds may produce too few candidates in some outer training splits.
The recommended response is to justify and relax thresholds. For an exploratory
run, `allow_fallback = TRUE` fills the pool with the best finite ranked markers;
every affected split is marked by `used_fallback`.

## Important interpretation

`result$nested$outer_summary` contains the internal validation result. The
training AUC of `result$final_model` is apparent performance and must not replace
nested or external validation.

The package does not claim that a selected panel is clinically validated. A
separate cohort, locked preprocessing procedure, and prospectively specified
threshold are required before clinical use.

## Before publishing to CRAN

Edit the placeholder maintainer identity, email, GitHub URL, and BugReports URL
in `DESCRIPTION`, then run:

```r
install.packages(c("devtools", "roxygen2", "testthat"))
devtools::document()
devtools::test()
devtools::check()
```

## Author

Fuhao Jiang <emr39515@gmail.com>
