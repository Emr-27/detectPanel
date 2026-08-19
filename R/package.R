#' detectPanel: leakage-aware biomarker panel discovery
#'
#' `detectPanel` discovers small biomarker panels for binary outcomes. Its
#' nested-validation workflow repeats marker preselection and panel search
#' inside each outer training split, which avoids evaluating a model on data
#' that were already used to choose its features.
#'
#' @author Fuhao Jiang <emr39515@gmail.com>
#' @name detectPanel
#' @docType package
NULL

utils::globalVariables(c(
  "change", "expression", "feature", "fpr", "frequency", "group",
  "log2FoldChange", "neg_log10_padj", "panel", "sensitivity"
))
