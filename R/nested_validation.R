.nested_stage <- function(stage, expr) {
  tryCatch(
    force(expr),
    error = function(e) {
      stop(paste0(stage, ": ", conditionMessage(e)), call. = FALSE)
    }
  )
}

.training_candidate_pool <- function(
    train_counts,
    train_outcome,
    transform_method,
    prior_count,
    detection_threshold,
    detectability_weights,
    preselection_weights,
    min_mean,
    min_median,
    min_detection,
    min_group_detection,
    min_auc,
    candidate_n,
    panel_size,
    allow_fallback
) {
  detectability <- score_detectability(
    train_counts,
    group = factor(train_outcome, levels = c(0, 1)),
    detection_threshold = detection_threshold,
    weights = detectability_weights
  )
  train_expression <- transform_assay(
    train_counts, method = transform_method, prior_count = prior_count
  )
  marker_metrics <- evaluate_markers(train_expression, train_outcome)
  selected <- preselect_markers(
    detectability = detectability,
    marker_metrics = marker_metrics,
    min_mean = min_mean,
    min_median = min_median,
    min_detection = min_detection,
    min_group_detection = min_group_detection,
    min_auc = min_auc,
    weights = preselection_weights,
    fallback_n = if (isTRUE(allow_fallback)) max(candidate_n, panel_size) else 0L
  )
  selected <- utils::head(selected, candidate_n)
  list(
    candidates = selected,
    detectability = detectability,
    marker_metrics = marker_metrics,
    expression = train_expression,
    audit = attr(selected, "audit")
  )
}

#' Nested cross-validation of biomarker panel discovery
#'
#' Candidate scoring, hard filtering, panel combination search, scaling, model
#' fitting, and threshold estimation are performed using only each outer
#' analysis set. Outer assessment samples are used only for final prediction.
#' If all outer splits fail, the raised error includes split-specific,
#' stage-specific diagnostic messages.
#'
#' @param counts Non-negative feature-by-sample count matrix.
#' @param metadata Sample metadata.
#' @param outcome Name of the binary outcome column.
#' @param positive Positive-class label.
#' @param sample_id Optional metadata sample-ID column.
#' @param exclude Exact feature names to remove before resampling.
#' @param exclude_regex Optional name pattern to remove before resampling.
#' @param panel_size Number of markers in each panel.
#' @param candidate_n Maximum number of training-selected candidates searched.
#' @param transform_method Assay transformation used for modeling.
#' @param prior_count Offset for log transformations.
#' @param detection_threshold Raw-count detection threshold.
#' @param detectability_weights Weights for abundance, detection, and stability.
#' @param preselection_weights Weights for detectability and AUC ranks.
#' @param min_mean,min_median,min_detection,min_group_detection,min_auc Hard
#'   training-only candidate filters.
#' @param allow_fallback If `TRUE`, fill an undersized candidate pool with the
#'   highest-ranked finite markers and flag the affected outer split.
#' @param outer_v,outer_repeats Outer cross-validation settings.
#' @param inner_v,inner_repeats Inner panel-search settings.
#' @param near_best_tolerance Predictive tolerance used for panel tie breaking.
#' @param seed Random seed.
#' @param max_combinations Safety limit per inner panel search.
#' @return A `detectPanel_nested` object.
#' @export
nested_validate_panels <- function(
    counts,
    metadata,
    outcome,
    positive = NULL,
    sample_id = NULL,
    exclude = character(),
    exclude_regex = NULL,
    panel_size = 3L,
    candidate_n = 10L,
    transform_method = c("log2_cpm", "log2", "none"),
    prior_count = 1,
    detection_threshold = 10,
    detectability_weights = c(0.50, 0.35, 0.15),
    preselection_weights = c(0.50, 0.50),
    min_mean = 100,
    min_median = 20,
    min_detection = 0.80,
    min_group_detection = 0.50,
    min_auc = 0.80,
    allow_fallback = FALSE,
    outer_v = 5L,
    outer_repeats = 2L,
    inner_v = 5L,
    inner_repeats = 3L,
    near_best_tolerance = 0.01,
    seed = 123L,
    max_combinations = 50000L
) {
  transform_method <- match.arg(transform_method)
  validated <- validate_assay_data(
    counts, metadata, outcome, positive = positive, sample_id = sample_id
  )
  excluded <- exclude_features(
    validated$assay, features = exclude, regex = exclude_regex
  )
  assay <- excluded$assay
  y <- validated$outcome
  sample_names <- colnames(assay)
  panel_size <- as.integer(panel_size)[1L]
  candidate_n <- as.integer(candidate_n)[1L]
  if (candidate_n < panel_size) {
    stop("`candidate_n` must be at least `panel_size`.", call. = FALSE)
  }
  outer_splits <- make_repeated_stratified_splits(
    y, v = outer_v, repeats = outer_repeats, seed = seed
  )

  prediction_parts <- list()
  metric_parts <- list()
  split_parts <- list()
  candidate_tables <- list()
  searches <- list()
  error_parts <- list()

  for (i in seq_along(outer_splits)) {
    split <- outer_splits[[i]]
    train <- split$analysis
    test <- split$assessment
    split_id <- sprintf("R%02dF%02d", split$repeat_id, split$fold)

    result <- tryCatch({
      pool <- .nested_stage("candidate selection", .training_candidate_pool(
        train_counts = assay[, train, drop = FALSE],
        train_outcome = y[train],
        transform_method = transform_method,
        prior_count = prior_count,
        detection_threshold = detection_threshold,
        detectability_weights = detectability_weights,
        preselection_weights = preselection_weights,
        min_mean = min_mean,
        min_median = min_median,
        min_detection = min_detection,
        min_group_detection = min_group_detection,
        min_auc = min_auc,
        candidate_n = candidate_n,
        panel_size = panel_size,
        allow_fallback = allow_fallback
      ))
      if (nrow(pool$candidates) < panel_size) {
        stop(
          "candidate selection: Training split produced too few finite candidate markers.",
          call. = FALSE
        )
      }
      inner_splits <- .nested_stage("inner resampling", make_repeated_stratified_splits(
        y[train], v = inner_v, repeats = inner_repeats,
        seed = seed + 1000L + i
      ))
      search <- .nested_stage("inner panel search", search_panels(
        expression = pool$expression,
        outcome = y[train],
        candidates = pool$candidates$feature,
        panel_size = panel_size,
        splits = inner_splits,
        candidate_table = pool$candidates,
        near_best_tolerance = near_best_tolerance,
        max_combinations = max_combinations
      ))
      model <- .nested_stage("outer-training model fit", fit_panel(
        assay = assay[, train, drop = FALSE],
        outcome = y[train],
        features = search$best_panel,
        transform_method = transform_method,
        prior_count = prior_count,
        threshold_method = "youden"
      ))
      probability <- .nested_stage("outer-assessment prediction", predict(
        model, assay[, test, drop = FALSE], type = "response"
      ))
      if (length(probability) != length(test) || any(!is.finite(probability))) {
        stop("outer-assessment prediction: non-finite or misaligned probabilities.",
             call. = FALSE)
      }
      predicted <- as.integer(probability >= model$threshold)
      metrics <- .classification_metrics(y[test], probability, model$threshold)
      list(pool = pool, search = search, model = model,
           probability = probability, predicted = predicted,
           metrics = metrics)
    }, error = function(e) e)

    if (inherits(result, "error")) {
      error_parts[[length(error_parts) + 1L]] <- data.frame(
        split_id = split_id,
        repeat_id = split$repeat_id,
        fold = split$fold,
        message = conditionMessage(result),
        stringsAsFactors = FALSE
      )
      next
    }

    panel <- result$search$best_panel
    prediction_parts[[length(prediction_parts) + 1L]] <- data.frame(
      sample = sample_names[test],
      observed = y[test],
      probability = as.numeric(result$probability),
      predicted = result$predicted,
      threshold = result$model$threshold,
      split_id = split_id,
      repeat_id = split$repeat_id,
      fold = split$fold,
      panel = .panel_key(panel),
      stringsAsFactors = FALSE
    )
    metric_parts[[length(metric_parts) + 1L]] <- cbind(
      data.frame(split_id = split_id, repeat_id = split$repeat_id,
                 fold = split$fold, panel = .panel_key(panel),
                 stringsAsFactors = FALSE),
      result$metrics
    )
    split_parts[[length(split_parts) + 1L]] <- data.frame(
      split_id = split_id,
      repeat_id = split$repeat_id,
      fold = split$fold,
      panel = .panel_key(panel),
      inner_mean_CV_AUC = result$search$near_best$mean_CV_AUC[1L],
      candidate_count = nrow(result$pool$candidates),
      used_fallback = any(result$pool$candidates$selected_by_fallback),
      stringsAsFactors = FALSE
    )
    candidate_tables[[split_id]] <- result$pool$candidates
    searches[[split_id]] <- result$search
  }

  predictions <- .bind_rows_base(prediction_parts)
  split_metrics <- .bind_rows_base(metric_parts)
  split_selections <- .bind_rows_base(split_parts)
  errors <- .bind_rows_base(error_parts)
  if (!nrow(predictions)) {
    details <- if (nrow(errors)) {
      paste(sprintf("  %s: %s", errors$split_id, errors$message),
            collapse = "\n")
    } else {
      "  No split-specific error messages were recorded."
    }
    stop(
      paste0(
        "Every outer split failed. Split-specific errors:\n",
        details
      ),
      call. = FALSE
    )
  }

  aggregated <- stats::aggregate(
    probability ~ sample + observed,
    data = predictions,
    FUN = mean
  )
  vote <- stats::aggregate(
    predicted ~ sample + observed,
    data = predictions,
    FUN = function(x) as.integer(mean(x) >= 0.5)
  )
  aggregated <- merge(aggregated, vote, by = c("sample", "observed"),
                      all.x = TRUE, sort = FALSE)
  pooled_auc <- .auc_fixed(aggregated$observed, aggregated$probability)

  panel_count <- stats::aggregate(
    split_id ~ panel, data = split_selections, length
  )
  names(panel_count)[2L] <- "count"
  panel_auc <- stats::aggregate(
    inner_mean_CV_AUC ~ panel, data = split_selections, mean
  )
  panel_frequency <- merge(panel_count, panel_auc, by = "panel", sort = FALSE)
  panel_frequency$frequency <- panel_frequency$count / nrow(split_selections)
  panel_frequency <- panel_frequency[
    order(-panel_frequency$count, -panel_frequency$inner_mean_CV_AUC,
          panel_frequency$panel),
    , drop = FALSE
  ]
  rownames(panel_frequency) <- NULL
  feature_vector <- unlist(strsplit(split_selections$panel, " \\+ "))
  feature_tab <- sort(table(feature_vector), decreasing = TRUE)
  feature_frequency <- data.frame(
    feature = names(feature_tab),
    count = as.integer(feature_tab),
    frequency = as.integer(feature_tab) / nrow(split_selections),
    stringsAsFactors = FALSE
  )
  auc_values <- split_metrics$AUC[is.finite(split_metrics$AUC)]
  outer_summary <- data.frame(
    requested_splits = length(outer_splits),
    valid_splits = nrow(split_metrics),
    failed_splits = nrow(errors),
    mean_outer_AUC = if (length(auc_values)) mean(auc_values) else NA_real_,
    sd_outer_AUC = if (length(auc_values) > 1L) stats::sd(auc_values) else 0,
    median_outer_AUC = if (length(auc_values)) stats::median(auc_values) else NA_real_,
    minimum_outer_AUC = if (length(auc_values)) min(auc_values) else NA_real_,
    maximum_outer_AUC = if (length(auc_values)) max(auc_values) else NA_real_,
    pooled_sample_level_AUC = pooled_auc
  )

  out <- list(
    predictions = predictions,
    aggregated_predictions = aggregated,
    split_metrics = split_metrics,
    split_selections = split_selections,
    panel_frequency = panel_frequency,
    feature_frequency = feature_frequency,
    outer_summary = outer_summary,
    candidate_tables = candidate_tables,
    searches = searches,
    errors = errors,
    excluded = excluded$excluded,
    positive = validated$positive,
    negative = validated$negative,
    settings = list(
      panel_size = panel_size,
      candidate_n = candidate_n,
      transform_method = transform_method,
      prior_count = prior_count,
      detection_threshold = detection_threshold,
      min_mean = min_mean,
      min_median = min_median,
      min_detection = min_detection,
      min_group_detection = min_group_detection,
      min_auc = min_auc,
      allow_fallback = allow_fallback,
      outer_v = attr(outer_splits, "v"),
      outer_repeats = outer_repeats,
      inner_v = inner_v,
      inner_repeats = inner_repeats,
      seed = seed
    )
  )
  class(out) <- "detectPanel_nested"
  out
}

#' @export
print.detectPanel_nested <- function(x, ...) {
  cat("detectPanel nested validation\n")
  cat("  Valid outer splits:", x$outer_summary$valid_splits,
      "of", x$outer_summary$requested_splits, "\n")
  cat("  Mean outer AUC:", sprintf("%.3f", x$outer_summary$mean_outer_AUC), "\n")
  cat("  Sample-level pooled AUC:",
      sprintf("%.3f", x$outer_summary$pooled_sample_level_AUC), "\n")
  if (nrow(x$panel_frequency)) {
    cat("  Most frequent panel:", x$panel_frequency$panel[1L], "\n")
  }
  invisible(x)
}

#' @export
summary.detectPanel_nested <- function(object, ...) {
  list(
    outer_summary = object$outer_summary,
    panel_frequency = object$panel_frequency,
    feature_frequency = object$feature_frequency,
    failed_splits = object$errors
  )
}

#' @export
plot.detectPanel_nested <- function(x, type = c("roc", "feature_frequency",
                                                 "panel_frequency"), ...) {
  type <- match.arg(type)
  if (type == "roc") {
    roc <- .roc_points(x$aggregated_predictions$observed,
                       x$aggregated_predictions$probability)
    return(
      ggplot2::ggplot(roc, ggplot2::aes(x = fpr, y = sensitivity)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        ggplot2::coord_equal() +
        ggplot2::labs(
          title = "Nested out-of-fold ROC",
          subtitle = sprintf("Sample-level pooled AUC = %.3f",
                             x$outer_summary$pooled_sample_level_AUC),
          x = "1 - Specificity", y = "Sensitivity"
        ) +
        ggplot2::theme_bw()
    )
  }
  if (type == "feature_frequency") {
    d <- x$feature_frequency
    d$feature <- factor(d$feature, levels = rev(d$feature))
    return(
      ggplot2::ggplot(d, ggplot2::aes(x = feature, y = frequency)) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::labs(title = "Feature selection frequency",
                      x = NULL, y = "Outer-split frequency") +
        ggplot2::theme_bw()
    )
  }
  d <- utils::head(x$panel_frequency, 15L)
  d$panel <- factor(d$panel, levels = rev(d$panel))
  ggplot2::ggplot(d, ggplot2::aes(x = panel, y = frequency)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Panel selection frequency",
                  x = NULL, y = "Outer-split frequency") +
    ggplot2::theme_bw()
}
