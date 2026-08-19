# detectPanel 0.1.0

## Initial release

`detectPanel` provides a leakage-aware workflow for discovery and
validation of small biomarker panels.

## Features

-   Detectability-aware feature filtering.
-   Small-panel search with nested cross-validation.
-   Logistic panel fitting and prediction.
-   L2-penalized fallback for unstable logistic regression.
-   Prediction support for fitted models.
-   Panel and feature stability summaries.
-   Export and visualization utilities.

## Validation

Tested using simulated datasets, unit tests, `R CMD check --as-cran`,
and an independent GEO miRNA-seq workflow based on GSE185062.
