\name{flag_pca_outliers}
\alias{flag_pca_outliers}
\title{Flag Potential PCA Outliers Without Deleting Samples}
\usage{
flag_pca_outliers(expression, group = NULL, n_components = 5L,
  threshold = NULL)
}
\arguments{
  \item{expression}{Transformed feature-by-sample matrix.}
  \item{group}{Optional sample group vector.}
  \item{n_components}{Maximum number of principal components used.}
  \item{threshold}{Optional positive robust-distance threshold.}
}
\value{A sample-level data frame with PCA coordinates, robust distance, and a flag.}
\description{Provides non-destructive PCA quality-control flags for review alongside independent QC evidence.}
