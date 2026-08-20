# detectPanel

`detectPanel` is an R package for discovering, validating, and applying small
biomarker panels from feature-by-sample count or expression matrices.

It provides an end-to-end workflow for feature filtering, small-panel search,
nested cross-validation, model fitting, prediction on new samples, stability
assessment, visualization, and export.

## Overview

High-dimensional molecular datasets often contain many candidate biomarkers
but relatively few samples. `detectPanel` is designed for binary-outcome
problems where the goal is to identify a compact panel while reducing
information leakage during model selection and validation.

A typical workflow is:

1. prepare a feature-by-sample assay matrix and sample metadata;
2. filter or exclude unsuitable features;
3. search for compact candidate panels;
4. estimate internal performance with nested cross-validation;
5. refit a final model on all available training samples;
6. apply the fitted model to new samples.

## Key features

- Detectability-aware feature filtering
- Feature exclusion by exact name or regular expression
- Small biomarker panel search
- Repeated stratified resampling
- Leakage-aware nested cross-validation
- Logistic panel fitting with stability safeguards
- Automatic L2-penalized fallback for unstable logistic regression fits
- Prediction from fitted `detectPanel_result` objects
- Panel and feature selection-frequency summaries
- Marker, ROC, heatmap, volcano, and PCA-QC plotting utilities
- Export of analysis results

## Installation

### From GitHub

```r
install.packages("remotes")
remotes::install_github("Emr-27/detectPanel")
```

### From CRAN

Once the package is available on CRAN:

```r
install.packages("detectPanel")
```

## Quick start

The example below generates a small synthetic count matrix and can be run
directly after installing the package.

```r
library(detectPanel)

set.seed(42)

n_samples <- 48
n_features <- 18

group <- rep(c("Control", "Case"), each = n_samples / 2)

counts <- matrix(
  rpois(n_features * n_samples, lambda = 80),
  nrow = n_features,
  dimnames = list(
    paste0("marker", seq_len(n_features)),
    paste0("S", seq_len(n_samples))
  )
)

# Add signal to a few features in the positive class
counts[1:3, group == "Case"] <-
  counts[1:3, group == "Case"] + 80

metadata <- data.frame(
  group = group,
  row.names = colnames(counts)
)

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
```

The selected panel is available from:

```r
result$final_panel
```

## Prediction

A fitted `detectPanel_result` object can be applied directly to a new
feature-by-sample matrix.

```r
new_counts <- counts[, 1:4, drop = FALSE]

probability <- predict(
  result,
  new_counts,
  type = "response"
)

classification <- predict(
  result,
  new_counts,
  type = "class"
)

probability
classification
```

New data must contain the features required by the fitted model. The saved
model retains the preprocessing and decision information needed for
prediction.

## Validation and interpretation

`detectPanel` uses nested cross-validation to estimate internal validation
performance while keeping candidate selection and panel search within the
training data of each outer split.

Nested validation summaries are available from:

```r
result$nested$outer_summary
```

The final model is refitted on all available training samples for downstream
prediction. Its training performance is apparent performance and should not be
used as a substitute for nested cross-validation or independent external
validation.

A selected panel should therefore be interpreted as a candidate biomarker
model rather than evidence of clinical validation. Independent cohorts and a
locked preprocessing and prediction procedure are recommended before
prospective or clinical use.

## Visualization

Depending on the fitted result and analysis stage, `detectPanel` provides
plotting functions for marker expression, marker ROC curves, panel heatmaps,
volcano plots, and PCA-based QC.

Examples include:

```r
plot(result, type = "roc")
plot(result, type = "feature_frequency")
```

See the function documentation for additional plotting options.

## Documentation

Detailed workflows are included as package vignettes:

```r
browseVignettes("detectPanel")
```

The package includes a vignette covering nested validation and the
recommended interpretation of internal validation results.

Function-level documentation is available through R help, for example:

```r
?discover_panel
?fit_panel
?nested_validate_panels
?export_results
```

## Citation

If you use `detectPanel` in research, please cite the package:

```r
citation("detectPanel")
```

## License

`detectPanel` is distributed under the MIT License.
