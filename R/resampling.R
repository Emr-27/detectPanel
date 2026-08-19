#' Create repeated stratified cross-validation splits
#'
#' The returned splits can be reused for every candidate panel, allowing fair
#' paired comparison between panels. Panel scoring requires a complete finite
#' out-of-fold prediction vector for each repeat; failed folds are not silently
#' omitted from AUC calculation.
#'
#' @param outcome Binary vector.
#' @param v Number of folds.
#' @param repeats Number of repeated fold assignments.
#' @param seed Random seed.
#' @return A list of split objects, each containing `analysis`, `assessment`,
#'   `repeat_id`, and `fold` indices.
#' @export
make_repeated_stratified_splits <- function(outcome, v = 5L, repeats = 10L,
                                             seed = 123L) {
  if (is.numeric(outcome) && all(outcome %in% c(0, 1))) {
    y <- factor(outcome, levels = c(0, 1))
  } else {
    y <- factor(outcome)
  }
  if (anyNA(y) || nlevels(y) != 2L) {
    stop("`outcome` must be a non-missing binary vector.", call. = FALSE)
  }
  v <- as.integer(v)[1L]
  repeats <- as.integer(repeats)[1L]
  if (!is.finite(v) || v < 2L || !is.finite(repeats) || repeats < 1L) {
    stop("`v` must be at least 2 and `repeats` at least 1.", call. = FALSE)
  }
  v_use <- min(v, min(table(y)))
  if (v_use < 2L) stop("Each class must contain at least two samples.",
                       call. = FALSE)
  set.seed(seed)
  splits <- list()
  counter <- 1L
  for (r in seq_len(repeats)) {
    fold_id <- integer(length(y))
    for (cls in levels(y)) {
      idx <- sample(which(y == cls))
      fold_id[idx] <- rep(seq_len(v_use), length.out = length(idx))
    }
    for (fold in seq_len(v_use)) {
      assessment <- which(fold_id == fold)
      analysis <- setdiff(seq_along(y), assessment)
      splits[[counter]] <- list(
        analysis = analysis,
        assessment = assessment,
        repeat_id = r,
        fold = fold
      )
      counter <- counter + 1L
    }
  }
  attr(splits, "v") <- v_use
  attr(splits, "repeats") <- repeats
  attr(splits, "seed") <- seed
  class(splits) <- c("detectPanel_splits", "list")
  splits
}

.evaluate_panel_on_splits <- function(expression, outcome, features, splits) {
  features <- as.character(features)
  if (!all(features %in% rownames(expression))) {
    return(data.frame(
      mean_CV_AUC = NA_real_, sd_CV_AUC = NA_real_,
      median_CV_AUC = NA_real_, minimum_CV_AUC = NA_real_,
      maximum_CV_AUC = NA_real_, valid_repeats = 0L,
      convergence_failures = length(splits)
    ))
  }
  repeat_ids <- sort(unique(vapply(splits, `[[`, integer(1), "repeat_id")))
  repeat_auc <- rep(NA_real_, length(repeat_ids))
  failures <- 0L

  for (r_i in seq_along(repeat_ids)) {
    current <- splits[vapply(splits, function(s) s$repeat_id == repeat_ids[r_i],
                             logical(1))]
    pred <- rep(NA_real_, ncol(expression))
    repeat_failed <- FALSE

    for (split in current) {
      train <- split$analysis
      test <- split$assessment
      x_train <- t(expression[features, train, drop = FALSE])
      x_test <- t(expression[features, test, drop = FALSE])
      scaling <- .scale_fit(x_train)
      x_test <- .scale_apply(x_test, scaling$center, scaling$scale)
      fitted <- .fit_logistic_matrix(scaling$x, outcome[train])
      if (is.null(fitted$fit) || !fitted$converged) {
        failures <- failures + 1L
        repeat_failed <- TRUE
        next
      }
      x_test_df <- as.data.frame(x_test, check.names = FALSE)
      names(x_test_df) <- fitted$predictor_names
      p <- tryCatch(
        as.numeric(stats::predict(fitted$fit, newdata = x_test_df,
                                  type = "response")),
        error = function(e) rep(NA_real_, length(test))
      )
      if (length(p) != length(test) || any(!is.finite(p))) {
        failures <- failures + 1L
        repeat_failed <- TRUE
        next
      }
      pred[test] <- p
    }

    # A repeated-CV AUC is valid only when every sample received a finite OOF
    # prediction. Never score a panel after silently dropping failed folds.
    if (!repeat_failed && all(is.finite(pred))) {
      repeat_auc[r_i] <- .auc_fixed(outcome, pred)
    }
  }

  valid <- repeat_auc[is.finite(repeat_auc)]
  if (!length(valid)) {
    return(data.frame(
      mean_CV_AUC = NA_real_, sd_CV_AUC = NA_real_,
      median_CV_AUC = NA_real_, minimum_CV_AUC = NA_real_,
      maximum_CV_AUC = NA_real_, valid_repeats = 0L,
      convergence_failures = failures
    ))
  }
  data.frame(
    mean_CV_AUC = mean(valid),
    sd_CV_AUC = if (length(valid) > 1L) stats::sd(valid) else 0,
    median_CV_AUC = stats::median(valid),
    minimum_CV_AUC = min(valid),
    maximum_CV_AUC = max(valid),
    valid_repeats = length(valid),
    convergence_failures = failures
  )
}
