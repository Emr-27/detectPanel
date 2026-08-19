# detectPanel 0.1.0

* Initial public release.
* Added detectability scoring based on abundance, detection rate, and robust
  within-group stability.
* Added shared repeated stratified folds for fair panel comparison.
* Added exhaustive small-panel search with logistic cross-validation.
* Added nested cross-validation that repeats candidate selection and panel
  search inside each outer training split.
* Added deployable fitted objects with saved transformation, centering,
  scaling, coefficients, and decision threshold.
* Added `predict()` support for fitted `detectPanel_result` objects.
* Added an internal L2-penalized logistic fallback for non-convergent or
  separation-prone ordinary logistic fits. Fitted objects record the actual
  fitting method.
* Tightened cross-validation scoring: incomplete out-of-fold predictions are
  not used for AUC calculation.
* Improved nested-validation failure diagnostics by reporting split-specific
  errors when all outer splits fail.
* Added regression tests covering prediction, separation handling, and
  nested-validation failure diagnostics.
