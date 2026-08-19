\name{validate_assay_data}
\alias{validate_assay_data}
\alias{exclude_features}
\alias{transform_assay}
\title{Validate, Filter, and Transform Assay Data}
\usage{
validate_assay_data(assay, metadata, outcome, positive = NULL, sample_id = NULL)
exclude_features(assay, features = character(), regex = NULL, ignore_case = TRUE)
transform_assay(assay, method = c("log2_cpm", "log2", "none"), prior_count = 1)
}
\arguments{
  \item{assay}{Numeric feature-by-sample matrix.}
  \item{metadata}{Sample metadata.}
  \item{outcome}{Binary outcome column name.}
  \item{positive}{Positive-class label.}
  \item{sample_id}{Optional metadata sample-ID column.}
  \item{features}{Exact feature names to remove.}
  \item{regex}{Optional feature-name regular expression.}
  \item{ignore_case}{Whether regular-expression matching ignores case.}
  \item{method}{Transformation method.}
  \item{prior_count}{Positive log offset.}
}
\value{Aligned data, a filtering audit, or a transformed matrix, depending on the function.}
\description{Core utilities for reproducible assay preparation without hard-coded dataset rules.}
