\name{plot_marker_expression}
\alias{plot_marker_expression}
\alias{plot_marker_roc}
\alias{plot_volcano}
\alias{plot_panel_heatmap}
\title{Diagnostic Plotting Helpers}
\usage{
plot_marker_expression(expression, metadata, features, outcome,
  sample_id = NULL)

plot_marker_roc(expression, outcome, features, positive = NULL)

plot_volcano(results, feature_col = "feature",
  lfc_col = "log2FoldChange", padj_col = "padj",
  padj_cut = 0.05, lfc_cut = 1)

plot_panel_heatmap(expression, features, metadata = NULL,
  annotation_cols = NULL, scale = "row", ...)
}
\arguments{
  \item{expression}{Transformed feature-by-sample matrix.}
  \item{metadata}{Sample metadata aligned by row names.}
  \item{features}{Features to display.}
  \item{outcome}{Outcome column name or a binary vector, depending on function.}
  \item{sample_id}{Optional metadata sample-ID column.}
  \item{positive}{Positive-class label.}
  \item{results}{Differential-expression result data frame.}
  \item{feature_col, lfc_col, padj_col}{Differential-expression column names.}
  \item{padj_cut, lfc_cut}{Volcano significance cutoffs.}
  \item{annotation_cols}{Metadata columns used as heatmap annotations.}
  \item{scale}{Heatmap scaling option.}
  \item{...}{Additional arguments passed to \code{pheatmap::pheatmap()}.}
}
\value{A ggplot object or the object returned by \code{pheatmap::pheatmap()}.}
\description{Creates expression, ROC, volcano, and optional heatmap views from package results.}
