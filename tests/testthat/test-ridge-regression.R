test_that("separation-prone logistic fits use a finite ridge fallback", {
  set.seed(101)
  y <- rep(0:1, each = 12)
  x <- cbind(
    signal = c(rep(-8, 12), rep(8, 12)),
    noise = stats::rnorm(24)
  )

  fitted <- detectPanel:::.fit_logistic_matrix(x, y)
  expect_true(fitted$converged)
  expect_identical(fitted$method, "ridge_fallback")
  expect_s3_class(fitted$fit, "detectPanel_ridge_fit")
  expect_true(all(is.finite(stats::coef(fitted$fit))))

  newdata <- as.data.frame(x, check.names = FALSE)
  names(newdata) <- fitted$predictor_names
  p <- stats::predict(fitted$fit, newdata = newdata, type = "response")
  expect_length(p, length(y))
  expect_true(all(is.finite(p)))
  expect_true(all(p >= 0 & p <= 1))
})

test_that("all-outer-split failures expose split-specific diagnostics", {
  set.seed(102)
  y <- rep(c("Control", "Case"), each = 8)
  counts <- matrix(
    stats::rpois(2 * length(y), lambda = 50),
    nrow = 2,
    dimnames = list(c("m1", "m2"), paste0("s", seq_along(y)))
  )
  metadata <- data.frame(group = y, row.names = colnames(counts))

  expect_error(
    nested_validate_panels(
      counts = counts,
      metadata = metadata,
      outcome = "group",
      positive = "Case",
      panel_size = 3,
      candidate_n = 3,
      detection_threshold = 1,
      min_mean = 0,
      min_median = 0,
      min_detection = 0,
      min_group_detection = 0,
      min_auc = 0.5,
      allow_fallback = TRUE,
      outer_v = 2,
      outer_repeats = 1,
      inner_v = 2,
      inner_repeats = 1,
      seed = 11
    ),
    regexp = "R01F01",
    fixed = FALSE
  )
})

test_that("strong-signal nested validation completes with finite predictions", {
  set.seed(103)
  n <- 48L
  y <- rep(c("Control", "Case"), each = n / 2L)
  counts <- matrix(
    stats::rnbinom(10 * n, mu = 100, size = 20),
    nrow = 10,
    dimnames = list(paste0("m", seq_len(10)), paste0("s", seq_len(n)))
  )
  counts[1, y == "Case"] <- counts[1, y == "Case"] + 450
  counts[2, y == "Case"] <- counts[2, y == "Case"] + 300
  counts[3, y == "Case"] <- pmax(0, counts[3, y == "Case"] - 80)
  metadata <- data.frame(group = y, row.names = colnames(counts))

  fit <- nested_validate_panels(
    counts = counts,
    metadata = metadata,
    outcome = "group",
    positive = "Case",
    panel_size = 3,
    candidate_n = 6,
    transform_method = "log2_cpm",
    detection_threshold = 5,
    min_mean = 5,
    min_median = 2,
    min_detection = 0.2,
    min_group_detection = 0.1,
    min_auc = 0.5,
    allow_fallback = TRUE,
    outer_v = 3,
    outer_repeats = 1,
    inner_v = 3,
    inner_repeats = 1,
    seed = 17
  )

  expect_s3_class(fit, "detectPanel_nested")
  expect_gt(nrow(fit$predictions), 0)
  expect_true(all(is.finite(fit$predictions$probability)))
  expect_true(is.finite(fit$outer_summary$pooled_sample_level_AUC))
})


test_that("incomplete out-of-fold repeats are not scored", {
  expression <- matrix(
    c(1, 2, 8, 9,
      2, 1, 9, 8),
    nrow = 2, byrow = TRUE,
    dimnames = list(c("m1", "m2"), paste0("s", 1:4))
  )
  y <- c(0L, 0L, 1L, 1L)
  bad_splits <- list(
    list(analysis = c(1L, 2L), assessment = c(3L, 4L),
         repeat_id = 1L, fold = 1L),
    list(analysis = c(3L, 4L), assessment = c(1L, 2L),
         repeat_id = 1L, fold = 2L)
  )

  cv <- detectPanel:::.evaluate_panel_on_splits(
    expression = expression,
    outcome = y,
    features = c("m1", "m2"),
    splits = bad_splits
  )

  expect_equal(cv$valid_repeats, 0L)
  expect_true(is.na(cv$mean_CV_AUC))
  expect_equal(cv$convergence_failures, 2L)
})
