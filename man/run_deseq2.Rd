\name{run_deseq2}
\alias{run_deseq2}
\title{Optional DESeq2 Differential-Expression Adapter}
\usage{
run_deseq2(counts, metadata, design, contrast, min_total_count = 10, ...)
}
\arguments{
  \item{counts}{Feature-by-sample count matrix.}
  \item{metadata}{Sample metadata aligned by row names.}
  \item{design}{Design formula.}
  \item{contrast}{Contrast passed to \code{DESeq2::results()}.}
  \item{min_total_count}{Minimum total feature count.}
  \item{...}{Additional arguments passed to \code{DESeq2::results()}.}
}
\value{A list containing a fitted DESeq2 object and result data frame.}
\description{Produces differential-expression evidence separately from predictive validation.}
