# Internal helpers ---------------------------------------------------------

.as_numeric_matrix <- function(x, arg = "x") {
  if (is.data.frame(x)) {
    x <- as.matrix(x)
  }
  if (!is.matrix(x) || !is.numeric(x)) {
    stop(sprintf("`%s` must be a numeric matrix or numeric data frame.", arg),
         call. = FALSE)
  }
  storage.mode(x) <- "double"
  x
}

.safe_percent_rank <- function(x) {
  x <- as.numeric(x)
  out <- rep(NA_real_, length(x))
  ok <- is.finite(x)
  n_ok <- sum(ok)
  if (n_ok == 0L) return(out)
  if (n_ok == 1L) {
    out[ok] <- 1
    return(out)
  }
  out[ok] <- (rank(x[ok], ties.method = "average") - 1) / (n_ok - 1)
  out
}

.safe_mad <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (length(x) < 2L) return(NA_real_)
  stats::mad(x, constant = 1, na.rm = TRUE)
}

.safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

.safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else min(x)
}

.check_binary_outcome <- function(y, positive = NULL) {
  if (anyNA(y)) stop("Outcome contains missing values.", call. = FALSE)
  lev <- unique(as.character(y))
  if (length(lev) != 2L) {
    stop("Outcome must contain exactly two classes.", call. = FALSE)
  }
  if (is.null(positive)) positive <- lev[2L]
  positive <- as.character(positive)[1L]
  if (!positive %in% lev) {
    stop("`positive` is not present in the outcome.", call. = FALSE)
  }
  negative <- setdiff(lev, positive)[1L]
  list(
    y = as.integer(as.character(y) == positive),
    positive = positive,
    negative = negative
  )
}

.auc_fixed <- function(response, predictor) {
  ok <- is.finite(response) & is.finite(predictor)
  response <- as.integer(response[ok])
  predictor <- as.numeric(predictor[ok])
  if (length(response) < 4L || length(unique(response)) != 2L ||
      length(unique(predictor)) < 2L) {
    return(NA_real_)
  }
  n1 <- sum(response == 1L)
  n0 <- sum(response == 0L)
  ranks <- rank(predictor, ties.method = "average")
  (sum(ranks[response == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

.roc_points <- function(response, predictor) {
  ok <- is.finite(response) & is.finite(predictor)
  response <- as.integer(response[ok])
  predictor <- as.numeric(predictor[ok])
  if (!length(response) || length(unique(response)) != 2L) {
    return(data.frame(fpr = numeric(), sensitivity = numeric(),
                      specificity = numeric(), threshold = numeric()))
  }
  thresholds <- c(Inf, sort(unique(predictor), decreasing = TRUE), -Inf)
  n_pos <- sum(response == 1L)
  n_neg <- sum(response == 0L)
  rows <- lapply(thresholds, function(th) {
    pred <- as.integer(predictor >= th)
    tp <- sum(pred == 1L & response == 1L)
    fp <- sum(pred == 1L & response == 0L)
    data.frame(
      fpr = fp / n_neg,
      sensitivity = tp / n_pos,
      specificity = 1 - fp / n_neg,
      threshold = th
    )
  })
  do.call(rbind, rows)
}

.best_youden_threshold <- function(response, predictor) {
  pts <- .roc_points(response, predictor)
  if (!nrow(pts)) return(0.5)
  score <- pts$sensitivity + pts$specificity - 1
  keep <- is.finite(pts$threshold) & is.finite(score)
  if (!any(keep)) return(0.5)
  pts$threshold[which.max(replace(score, !keep, -Inf))]
}

.classification_metrics <- function(response, probability, threshold = 0.5) {
  ok <- is.finite(response) & is.finite(probability)
  response <- as.integer(response[ok])
  probability <- as.numeric(probability[ok])
  if (!length(response)) {
    return(data.frame(
      n = 0L, AUC = NA_real_, threshold = threshold,
      sensitivity = NA_real_, specificity = NA_real_, accuracy = NA_real_,
      balanced_accuracy = NA_real_, TN = NA_integer_, FP = NA_integer_,
      FN = NA_integer_, TP = NA_integer_
    ))
  }
  pred <- as.integer(probability >= threshold)
  tn <- sum(response == 0L & pred == 0L)
  fp <- sum(response == 0L & pred == 1L)
  fn <- sum(response == 1L & pred == 0L)
  tp <- sum(response == 1L & pred == 1L)
  sens <- if ((tp + fn) > 0L) tp / (tp + fn) else NA_real_
  spec <- if ((tn + fp) > 0L) tn / (tn + fp) else NA_real_
  data.frame(
    n = length(response),
    AUC = .auc_fixed(response, probability),
    threshold = threshold,
    sensitivity = sens,
    specificity = spec,
    accuracy = (tp + tn) / length(response),
    balanced_accuracy = mean(c(sens, spec), na.rm = TRUE),
    TN = tn, FP = fp, FN = fn, TP = tp
  )
}

.ridge_logistic_loss <- function(beta, design, y, lambda) {
  eta <- as.numeric(design %*% beta)
  softplus <- pmax(eta, 0) + log1p(exp(-abs(eta)))
  penalty <- 0.5 * lambda * sum(beta[-1L]^2)
  sum(softplus - y * eta) + penalty
}

.ridge_logistic_fit <- function(x, y, lambda = 1, maxit = 100L,
                                tolerance = 1e-8) {
  x <- as.matrix(x)
  y <- as.integer(y)
  predictor_names <- colnames(x)
  if (is.null(predictor_names)) {
    predictor_names <- paste0("x", seq_len(ncol(x)))
    colnames(x) <- predictor_names
  }
  design <- cbind(`(Intercept)` = 1, x)
  beta <- rep(0, ncol(design))
  prevalence <- mean(y)
  if (is.finite(prevalence) && prevalence > 0 && prevalence < 1) {
    beta[1L] <- stats::qlogis(prevalence)
  }
  lambda <- as.numeric(lambda)[1L]
  if (!is.finite(lambda) || lambda <= 0) lambda <- 1
  maxit <- max(1L, as.integer(maxit)[1L])
  converged <- FALSE
  iter_used <- 0L

  for (iter in seq_len(maxit)) {
    iter_used <- iter
    eta <- as.numeric(design %*% beta)
    prob <- stats::plogis(eta)
    weight <- pmax(prob * (1 - prob), 1e-8)
    gradient <- as.numeric(crossprod(design, prob - y))
    gradient[-1L] <- gradient[-1L] + lambda * beta[-1L]
    hessian <- crossprod(design, design * weight)
    hessian <- hessian + diag(c(0, rep(lambda, ncol(design) - 1L)))

    step <- tryCatch(
      as.numeric(solve(hessian, gradient)),
      error = function(e) {
        jitter <- diag(1e-7, nrow(hessian))
        tryCatch(as.numeric(solve(hessian + jitter, gradient)),
                 error = function(e2) rep(NA_real_, length(beta)))
      }
    )
    if (any(!is.finite(step))) break

    old_loss <- .ridge_logistic_loss(beta, design, y, lambda)
    scale_step <- 1
    accepted <- FALSE
    beta_new <- beta
    for (ls in seq_len(25L)) {
      proposal <- beta - scale_step * step
      new_loss <- .ridge_logistic_loss(proposal, design, y, lambda)
      if (is.finite(new_loss) && new_loss <= old_loss + 1e-12) {
        beta_new <- proposal
        accepted <- TRUE
        break
      }
      scale_step <- scale_step / 2
    }
    if (!accepted) break

    delta <- max(abs(beta_new - beta))
    beta <- beta_new
    if (is.finite(delta) && delta < tolerance) {
      converged <- TRUE
      break
    }
  }

  if (!converged && all(is.finite(beta))) {
    eta <- as.numeric(design %*% beta)
    prob <- stats::plogis(eta)
    gradient <- as.numeric(crossprod(design, prob - y))
    gradient[-1L] <- gradient[-1L] + lambda * beta[-1L]
    converged <- max(abs(gradient)) < 1e-5
  }

  coefficients <- stats::setNames(beta, c("(Intercept)", predictor_names))
  fit <- list(
    coefficients = coefficients,
    predictor_names = predictor_names,
    lambda = lambda,
    converged = isTRUE(converged),
    iterations = iter_used
  )
  class(fit) <- "detectPanel_ridge_fit"
  fit
}

#' @export
coef.detectPanel_ridge_fit <- function(object, ...) object$coefficients

#' @export
predict.detectPanel_ridge_fit <- function(object, newdata = NULL,
                                           type = c("link", "response"), ...) {
  type <- match.arg(type)
  if (is.null(newdata)) {
    stop("`newdata` is required for internal ridge predictions.", call. = FALSE)
  }
  x <- as.data.frame(newdata, check.names = FALSE)
  if (!all(object$predictor_names %in% names(x))) {
    stop("Ridge prediction data are missing required predictors.", call. = FALSE)
  }
  x <- as.matrix(x[, object$predictor_names, drop = FALSE])
  eta <- as.numeric(object$coefficients[1L] +
                    x %*% object$coefficients[-1L])
  if (type == "link") eta else stats::plogis(eta)
}

.fit_logistic_matrix <- function(x, y, ridge_lambda = 1) {
  x <- as.data.frame(x, check.names = FALSE)
  predictor_names <- paste0("x", seq_len(ncol(x)))
  names(x) <- predictor_names
  y <- as.integer(y)
  if (length(y) != nrow(x) || anyNA(y) || !all(y %in% c(0L, 1L)) ||
      length(unique(y)) != 2L) {
    return(list(
      fit = NULL,
      converged = FALSE,
      warnings = "Logistic training outcome must contain both 0 and 1 classes.",
      predictor_names = predictor_names,
      method = "failed",
      ridge_lambda = NA_real_
    ))
  }
  x$.outcome <- y
  warnings <- character()
  fit <- tryCatch(
    withCallingHandlers(
      stats::glm(.outcome ~ ., data = x, family = stats::binomial()),
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )

  glm_ok <- !inherits(fit, "error") && isTRUE(fit$converged) &&
    all(is.finite(stats::coef(fit)))
  separation_warning <- any(grepl(
    "did not converge|fitted probabilities numerically 0 or 1 occurred",
    warnings,
    ignore.case = TRUE
  ))
  extreme_glm <- FALSE
  if (glm_ok) {
    glm_coef <- stats::coef(fit)
    glm_prob <- stats::fitted(fit)
    extreme_glm <- any(abs(glm_coef[-1L]) > 25, na.rm = TRUE) ||
      any(glm_prob < 1e-8 | glm_prob > 1 - 1e-8, na.rm = TRUE)
  }
  separation_like <- separation_warning || extreme_glm

  if (glm_ok && !separation_like) {
    return(list(
      fit = fit,
      converged = TRUE,
      warnings = unique(warnings),
      predictor_names = predictor_names,
      method = "glm",
      ridge_lambda = NA_real_
    ))
  }

  glm_problem <- if (inherits(fit, "error")) {
    conditionMessage(fit)
  } else if (!isTRUE(fit$converged)) {
    "ordinary logistic regression did not converge"
  } else if (!all(is.finite(stats::coef(fit)))) {
    "ordinary logistic regression produced non-finite coefficients"
  } else {
    "ordinary logistic regression showed separation-like numerical behavior"
  }

  ridge_x <- x[, predictor_names, drop = FALSE]
  ridge_fit <- .ridge_logistic_fit(
    ridge_x,
    y = y,
    lambda = ridge_lambda
  )
  ridge_ok <- isTRUE(ridge_fit$converged) &&
    all(is.finite(stats::coef(ridge_fit)))
  fallback_message <- paste0(
    "Ordinary logistic fit was unstable (", glm_problem,
    "); used internal L2-penalized logistic fallback (lambda = ",
    format(ridge_fit$lambda, trim = TRUE), ")."
  )

  list(
    fit = if (ridge_ok) ridge_fit else NULL,
    converged = ridge_ok,
    warnings = unique(c(warnings, fallback_message)),
    predictor_names = predictor_names,
    method = if (ridge_ok) "ridge_fallback" else "failed",
    ridge_lambda = ridge_fit$lambda
  )
}

.scale_fit <- function(x) {
  x <- as.matrix(x)
  center <- colMeans(x, na.rm = TRUE)
  scale <- apply(x, 2L, stats::sd, na.rm = TRUE)
  scale[!is.finite(scale) | scale == 0] <- 1
  z <- sweep(sweep(x, 2L, center, "-"), 2L, scale, "/")
  list(x = z, center = center, scale = scale)
}

.scale_apply <- function(x, center, scale) {
  x <- as.matrix(x)
  sweep(sweep(x, 2L, center, "-"), 2L, scale, "/")
}

.validate_weights <- function(weights, n, arg) {
  weights <- as.numeric(weights)
  if (length(weights) != n || any(!is.finite(weights)) || any(weights < 0) ||
      abs(sum(weights) - 1) > 1e-8) {
    stop(sprintf("`%s` must contain %d non-negative values summing to 1.",
                 arg, n), call. = FALSE)
  }
  weights
}

.panel_key <- function(features) paste(sort(features), collapse = " + ")

.bind_rows_base <- function(parts) {
  parts <- parts[!vapply(parts, is.null, logical(1))]
  if (!length(parts)) return(data.frame())
  out <- do.call(rbind, parts)
  rownames(out) <- NULL
  out
}

#' Validate and align assay data
#'
#' @param assay Numeric feature-by-sample matrix.
#' @param metadata Data frame containing sample information.
#' @param outcome Name of the binary outcome column.
#' @param positive Positive-class label. Defaults to the second observed class.
#' @param sample_id Optional metadata column containing sample identifiers. If
#'   `NULL`, metadata row names are used.
#' @return A list containing aligned `assay`, `metadata`, binary `outcome`, and
#'   class labels.
#' @export
validate_assay_data <- function(assay, metadata, outcome, positive = NULL,
                                sample_id = NULL) {
  assay <- .as_numeric_matrix(assay, "assay")
  if (is.null(rownames(assay)) || anyDuplicated(rownames(assay))) {
    stop("Assay rows must have unique feature names.", call. = FALSE)
  }
  if (is.null(colnames(assay)) || anyDuplicated(colnames(assay))) {
    stop("Assay columns must have unique sample names.", call. = FALSE)
  }
  metadata <- as.data.frame(metadata, stringsAsFactors = FALSE)
  if (!is.null(sample_id)) {
    if (!sample_id %in% names(metadata)) {
      stop("`sample_id` column was not found in metadata.", call. = FALSE)
    }
    rownames(metadata) <- as.character(metadata[[sample_id]])
  }
  if (is.null(rownames(metadata)) || anyDuplicated(rownames(metadata))) {
    stop("Metadata must have unique sample row names or a `sample_id` column.",
         call. = FALSE)
  }
  if (!outcome %in% names(metadata)) {
    stop("`outcome` column was not found in metadata.", call. = FALSE)
  }
  common <- intersect(colnames(assay), rownames(metadata))
  if (length(common) < 4L) {
    stop("Fewer than four samples are shared by assay and metadata.",
         call. = FALSE)
  }
  assay <- assay[, common, drop = FALSE]
  metadata <- metadata[common, , drop = FALSE]
  cls <- .check_binary_outcome(metadata[[outcome]], positive)
  if (min(table(cls$y)) < 2L) {
    stop("Each class must contain at least two samples.", call. = FALSE)
  }
  list(
    assay = assay,
    metadata = metadata,
    outcome = cls$y,
    positive = cls$positive,
    negative = cls$negative,
    outcome_name = outcome
  )
}
