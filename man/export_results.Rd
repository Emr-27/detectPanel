\name{export_results}
\alias{export_results}
\title{Export detectPanel Results}
\usage{
export_results(x, path, save_model = TRUE)
}
\arguments{
  \item{x}{A \code{detectPanel_result} or \code{detectPanel_nested} object.}
  \item{path}{Output directory.}
  \item{save_model}{Whether to save an RDS object.}
}
\value{The normalized path, invisibly.}
\description{Writes nested metrics, predictions, selection frequencies, audits, and optional final-model files.}
