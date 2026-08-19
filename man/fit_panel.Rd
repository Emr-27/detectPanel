\name{fit_panel}
\alias{fit_panel}
\alias{predict.detectPanel_fit}
\alias{print.detectPanel_fit}
\alias{summary.detectPanel_fit}
\title{Fit and Predict with a Biomarker Panel}
\usage{
fit_panel(assay, outcome, features, positive = NULL,
  transform_method = c("log2_cpm", "log2", "none"), prior_count = 1,
  threshold_method = c("youden", "fixed"), fixed_threshold = 0.5)

\method{predict}{detectPanel_fit}(object, newdata,
  type = c("response", "class", "link"), ...)
}
\arguments{
  \item{assay}{Feature-by-sample count or expression matrix.}
  \item{outcome}{Binary outcome aligned to columns.}
  \item{features}{Features included in the logistic model.}
  \item{positive}{Positive-class label.}
  \item{transform_method}{Input transformation.}
  \item{prior_count}{Log offset.}
  \item{threshold_method}{Training threshold method.}
  \item{fixed_threshold}{Fixed probability threshold.}
  \item{object}{A fitted \code{detectPanel_fit} object.}
  \item{newdata}{Feature-by-new-sample matrix.}
  \item{type}{Prediction scale.}
  \item{...}{Additional arguments, currently unused.}
}
\value{A fitted model object or prediction vector. The fitted object records
\code{fit_method} and \code{ridge_lambda} when the internal ridge fallback is used.}
\description{Stores preprocessing and scaling parameters needed to predict new samples.
Ordinary logistic regression is attempted first. If it is non-convergent or
shows separation-like numerical behavior, an internal L2-penalized logistic
fit is used as a deterministic fallback.}
