#' Predict new samples from a fitted detectPanel result
#'
#' Applies the final model stored in a detectPanel result object.
#'
#' @param object A `detectPanel_result` object returned by [discover_panel()].
#' @param newdata A feature-by-sample numeric matrix.
#' @param type Prediction type: `"response"` for probabilities or `"class"`.
#' @param ... Additional arguments passed to the fitted model predictor.
#'
#' @return A numeric vector of probabilities or a factor of predicted classes.
#' @export
predict.detectPanel_result <- function(
    object,
    newdata,
    type = c("response", "class"),
    ...
) {
  if (!inherits(object, "detectPanel_result")) {
    stop("object must be a detectPanel_result object.",
         call. = FALSE)
  }

  if (is.null(object$final_model)) {
    stop("The supplied result object has no final_model.",
         call. = FALSE)
  }

  type <- match.arg(type)

  predict(
    object$final_model,
    newdata = newdata,
    type = type,
    ...
  )
}
