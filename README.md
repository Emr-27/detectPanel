# detectPanel

`detectPanel` is an R package for discovering and validating small
biomarker panels from count or expression matrices.

## Overview

The package provides marker filtering, small-panel discovery, nested
cross-validation, model fitting, and prediction on new samples.

## Features

-   Detectability-aware feature filtering
-   Small biomarker panel search
-   Nested cross-validation
-   Logistic panel fitting with stability safeguards
-   L2-penalized fallback for unstable logistic regression
-   Prediction from fitted models
-   Panel and feature stability summaries
-   Result export and visualization

## Installation

Development version:

``` r
install.packages("remotes")
remotes::install_github("Emr-27/detectPanel")
```

After CRAN release:

``` r
install.packages("detectPanel")
```

## Quick start

``` r
library(detectPanel)

result <- discover_panel(
  counts = counts,
  metadata = metadata,
  outcome = "group",
  positive = "Case",
  panel_size = 3,
  candidate_n = 8,
  seed = 42
)

predict(result, new_counts, type = "response")
```

## Validation

Nested cross-validation results are available from:

``` r
result$nested$outer_summary
```

Training performance should not replace nested validation or external
validation.

## Citation

``` r
citation("detectPanel")
```

## License

GPL (\>= 3)
