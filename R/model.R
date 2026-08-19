#' Fit a deployable logistic biomarker panel
#'
#' Ordinary logistic regression is attempted first. If that fit is
#' non-convergent or shows separation-like numerical behavior, detectPanel
#' uses a deterministic internal L2-penalized logistic fallback and records
#' the fitting method in the returned object.
#'
#' @param assay Feature-by-sample count or expression matrix.
#' @param outcome Binary vector aligned to assay columns.
#' @param features Features included in the model.
#' @param positive Positive-class label when outcome is not 0/1.
#' @param transform_method Transformation applied before fitting.
#' @param prior_count Offset used by [transform_assay()].
#' @param threshold_method `"youden"` estimates a threshold from training
#'   predictions; `"fixed"` uses `fixed_threshold`.
#' @param fixed_threshold Probability threshold used for fixed classification.
#' @return A `detectPanel_fit` object containing `fit_method` and, when used,
#'   `ridge_lambda`.
#' @export
fit_panel <- function(
    assay,
    outcome,
    features,
    positive = NULL,
    transform_method = c("log2_cpm", "log2", "none"),
    prior_count = 1,
    threshold_method = c("youden", "fixed"),
    fixed_threshold = 0.5
) {
  assay <- .as_numeric_matrix(assay, "assay")
  transform_method <- match.arg(transform_method)
  threshold_method <- match.arg(threshold_method)
  features <- unique(as.character(features))
  if (!length(features) || !all(features %in% rownames(assay))) {
    stop("All requested `features` must be present in assay row names.",
         call. = FALSE)
  }
  if (length(outcome) != ncol(assay)) {
    stop("`outcome` must be aligned to assay columns.", call. = FALSE)
  }
  if (is.numeric(outcome) && all(outcome %in% c(0, 1))) {
    y <- as.integer(outcome)
    positive_label <- "1"
    negative_label <- "0"
  } else {
    cls <- .check_binary_outcome(outcome, positive)
    y <- cls$y
    positive_label <- cls$positive
    negative_label <- cls$negative
  }
  expression <- transform_assay(assay, method = transform_method,
                                prior_count = prior_count)
  x <- t(expression[features, , drop = FALSE])
  scaling <- .scale_fit(x)
  fitted <- .fit_logistic_matrix(scaling$x, y)
  if (is.null(fitted$fit) || !fitted$converged) {
    stop(paste0(
      "Logistic model fitting failed after ordinary and internal ridge ",
      "attempts. Consider fewer features or a larger training set."
    ), call. = FALSE)
  }
  train_x <- as.data.frame(scaling$x, check.names = FALSE)
  names(train_x) <- fitted$predictor_names
  train_prob <- as.numeric(stats::predict(fitted$fit, newdata = train_x,
                                          type = "response"))
  threshold <- if (threshold_method == "youden") {
    .best_youden_threshold(y, train_prob)
  } else {
    as.numeric(fixed_threshold)[1L]
  }
  if (!is.finite(threshold) || threshold < 0 || threshold > 1) {
    stop("The classification threshold must be between 0 and 1.",
         call. = FALSE)
  }
  metrics <- .classification_metrics(y, train_prob, threshold)
  coefficients <- stats::coef(fitted$fit)
  coefficient_table <- data.frame(
    term = names(coefficients),
    feature = c("(Intercept)", features),
    coefficient = as.numeric(coefficients),
    stringsAsFactors = FALSE
  )
  out <- list(
    features = features,
    transform_method = transform_method,
    prior_count = prior_count,
    center = scaling$center,
    scale = scaling$scale,
    model = fitted$fit,
    predictor_names = fitted$predictor_names,
    fit_method = fitted$method,
    ridge_lambda = fitted$ridge_lambda,
    coefficient_table = coefficient_table,
    threshold = threshold,
    threshold_method = threshold_method,
    positive = positive_label,
    negative = negative_label,
    training_probability = train_prob,
    training_outcome = y,
    training_sample_names = colnames(assay),
    training_metrics = metrics,
    warnings = fitted$warnings
  )
  class(out) <- "detectPanel_fit"
  out
}

#' @export
predict.detectPanel_fit <- function(object, newdata, type = c("response", "class",
                                                               "link"), ...) {
  type <- match.arg(type)
  newdata <- .as_numeric_matrix(newdata, "newdata")
  if (is.null(rownames(newdata)) || !all(object$features %in% rownames(newdata))) {
    stop("`newdata` must contain all model features as row names.",
         call. = FALSE)
  }
  expression <- transform_assay(
    newdata,
    method = object$transform_method,
    prior_count = object$prior_count
  )
  x <- t(expression[object$features, , drop = FALSE])
  x <- .scale_apply(x, object$center, object$scale)
  x_df <- as.data.frame(x, check.names = FALSE)
  names(x_df) <- object$predictor_names
  if (type == "link") {
    return(as.numeric(stats::predict(object$model, newdata = x_df,
                                     type = "link")))
  }
  probability <- as.numeric(stats::predict(object$model, newdata = x_df,
                                            type = "response"))
  names(probability) <- colnames(newdata)
  if (type == "response") return(probability)
  factor(
    ifelse(probability >= object$threshold, object$positive, object$negative),
    levels = c(object$negative, object$positive)
  )
}

#' @export
print.detectPanel_fit <- function(x, ...) {
  cat("detectPanel fitted logistic panel\n")
  cat("  Features:", paste(x$features, collapse = " + "), "\n")
  cat("  Transformation:", x$transform_method, "\n")
  cat("  Fit method:", x$fit_method, "\n")
  cat("  Training AUC:", sprintf("%.3f", x$training_metrics$AUC), "\n")
  cat("  Threshold:", sprintf("%.3f", x$threshold), "\n")
  invisible(x)
}

#' @export
summary.detectPanel_fit <- function(object, ...) {
  list(
    features = object$features,
    coefficients = object$coefficient_table,
    training_metrics = object$training_metrics,
    threshold = object$threshold,
    fit_method = object$fit_method,
    ridge_lambda = object$ridge_lambda,
    warnings = object$warnings
  )
}
