#' Optional DESeq2 differential-expression adapter
#'
#' Differential expression is deliberately kept separate from the predictive
#' nested-validation workflow. It can be used as biological support, but it is
#' not required as a hard candidate filter.
#'
#' @param counts Non-negative feature-by-sample count matrix.
#' @param metadata Sample metadata aligned by row names.
#' @param design A model formula passed to `DESeq2::DESeqDataSetFromMatrix()`.
#' @param contrast Contrast passed to `DESeq2::results()`.
#' @param min_total_count Minimum total count retained before fitting.
#' @param ... Additional arguments passed to `DESeq2::results()`.
#' @return A list containing the fitted DESeq2 data set and a result data frame.
#' @export
run_deseq2 <- function(counts, metadata, design, contrast,
                       min_total_count = 10, ...) {
  if (!requireNamespace("DESeq2", quietly = TRUE)) {
    stop("Package 'DESeq2' is required for `run_deseq2()`.", call. = FALSE)
  }
  counts <- .as_numeric_matrix(counts, "counts")
  metadata <- as.data.frame(metadata)
  if (is.null(colnames(counts)) || is.null(rownames(metadata))) {
    stop("Counts need sample column names and metadata needs sample row names.",
         call. = FALSE)
  }
  common <- intersect(colnames(counts), rownames(metadata))
  if (!length(common)) stop("Counts and metadata have no shared samples.",
                            call. = FALSE)
  counts <- round(counts[, common, drop = FALSE])
  metadata <- metadata[common, , drop = FALSE]
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = counts,
    colData = metadata,
    design = design
  )
  keep <- rowSums(DESeq2::counts(dds)) >= min_total_count
  dds <- dds[keep, ]
  dds <- DESeq2::DESeq(dds)
  result <- DESeq2::results(dds, contrast = contrast, ...)
  result_df <- as.data.frame(result)
  result_df$feature <- rownames(result_df)
  result_df <- result_df[, c("feature", setdiff(names(result_df), "feature")),
                         drop = FALSE]
  rownames(result_df) <- NULL
  list(dds = dds, results = result_df)
}
