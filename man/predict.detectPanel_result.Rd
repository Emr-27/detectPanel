\name{predict.detectPanel_result}
\alias{predict.detectPanel_result}
\title{Predict samples using a fitted detectPanel result}
\usage{
\method{predict}{detectPanel_result}(object, newdata,
  type = c("response", "class"), ...)
}
\description{
Applies the final fitted panel model from a \code{detectPanel_result}
object to new count data.
}
\arguments{
\item{object}{A \code{detectPanel_result} object returned by
\code{\link{discover_panel}}.}
\item{newdata}{A feature-by-sample numeric matrix.}
\item{type}{Prediction output type: probabilities or classes.}
\item{...}{Additional arguments.}
}
\value{
A vector of predicted probabilities or classes.
}
\examples{
# See discover_panel examples.
}
\keyword{prediction}
