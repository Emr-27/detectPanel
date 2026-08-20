# detectPanel 0.1.0

## Initial release

`detectPanel` provides a leakage-aware workflow for discovery and validation
of small biomarker panels.

## Features

- Detectability-aware feature filtering.
- Small-panel search with nested cross-validation.
- Repeated stratified resampling utilities.
- Logistic panel fitting with automatic L2-penalized fallback for unstable
  or separation-prone fits.
- Prediction support for fitted `detectPanel_result` objects.
- Panel and feature selection-frequency summaries.
- Export and visualization utilities.

## Validation

The package is covered by simulated-data tests and unit tests, including
nested-validation, model-fitting, prediction, resampling, and plotting
workflows. Release candidates are checked with `R CMD check --as-cran`.
