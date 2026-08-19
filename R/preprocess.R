#' Exclude assay features by name or regular expression
#'
#' @param assay Numeric feature-by-sample matrix.
#' @param features Exact feature names to remove.
#' @param regex Optional regular expression matched against feature names.
#' @param ignore_case Whether regular-expression matching ignores case.
#' @return A list with filtered `assay` and an `excluded` audit table.
#' @export
exclude_features <- function(assay, features = character(), regex = NULL,
                             ignore_case = TRUE) {
  assay <- .as_numeric_matrix(assay, "assay")
  if (is.null(rownames(assay))) {
    stop("Assay must contain feature row names.", call. = FALSE)
  }
  exact <- rownames(assay) %in% features
  pattern <- rep(FALSE, nrow(assay))
  if (!is.null(regex) && nzchar(regex)) {
    pattern <- grepl(regex, rownames(assay), ignore.case = ignore_case,
                     perl = TRUE)
  }
  remove <- exact | pattern
  reason <- ifelse(exact & pattern, "exact_and_regex",
                   ifelse(exact, "exact", ifelse(pattern, "regex", NA_character_)))
  audit <- data.frame(
    feature = rownames(assay)[remove],
    reason = reason[remove],
    stringsAsFactors = FALSE
  )
  list(assay = assay[!remove, , drop = FALSE], excluded = audit)
}

#' Transform count or expression data
#'
#' @param assay Numeric feature-by-sample matrix.
#' @param method One of `"log2_cpm"`, `"log2"`, or `"none"`.
#' @param prior_count Positive offset used before log transformation.
#' @return Transformed feature-by-sample matrix.
#' @export
transform_assay <- function(assay, method = c("log2_cpm", "log2", "none"),
                            prior_count = 1) {
  assay <- .as_numeric_matrix(assay, "assay")
  method <- match.arg(method)
  if (!is.numeric(prior_count) || length(prior_count) != 1L ||
      !is.finite(prior_count) || prior_count < 0) {
    stop("`prior_count` must be a finite non-negative number.", call. = FALSE)
  }
  if (method == "none") return(assay)
  if (prior_count <= 0) {
    stop("`prior_count` must be positive for log transformations.",
         call. = FALSE)
  }
  if (any(assay < 0, na.rm = TRUE)) {
    stop("Log count transformations require non-negative values.", call. = FALSE)
  }
  if (method == "log2") return(log2(assay + prior_count))
  lib_size <- colSums(assay, na.rm = TRUE)
  if (any(!is.finite(lib_size) | lib_size <= 0)) {
    stop("Every sample must have a positive finite library size for CPM.",
         call. = FALSE)
  }
  cpm <- sweep(assay, 2L, lib_size / 1e6, "/")
  log2(cpm + prior_count)
}
