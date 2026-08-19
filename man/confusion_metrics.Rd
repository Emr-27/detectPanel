\name{confusion_metrics}
\alias{confusion_metrics}
\title{Compute Confusion-Matrix Metrics}
\usage{
confusion_metrics(response, probability, threshold = 0.5)
}
\arguments{
  \item{response}{Non-missing binary 0/1 outcome.}
  \item{probability}{Predicted probability for class 1.}
  \item{threshold}{Probability classification threshold.}
}
\value{A one-row data frame with AUC, sensitivity, specificity, accuracy, balanced accuracy, and confusion counts.}
\description{Computes fixed-threshold classification metrics without optimizing on the evaluation data.}
