test_that("fitted panel predicts new samples", {
  d <- make_synthetic_panel_data()
  fit <- fit_panel(
    d$counts, d$y, features = c("m1", "m2", "m3"), positive = "Case"
  )
  expect_s3_class(fit, "detectPanel_fit")
  p <- predict(fit, d$counts[, 1:5, drop = FALSE], type = "response")
  cls <- predict(fit, d$counts[, 1:5, drop = FALSE], type = "class")
  expect_length(p, 5)
  expect_true(all(p >= 0 & p <= 1))
  expect_length(cls, 5)
  expect_equal(fit$features, c("m1", "m2", "m3"))
})


test_that("fit_panel records the fitting method and finite predictions", {
  d <- make_synthetic_panel_data(n = 36, p = 8, seed = 44)
  fit <- fit_panel(
    d$counts, d$y, features = rownames(d$counts)[1:3],
    positive = "Case", transform_method = "log2_cpm"
  )
  expect_true(fit$fit_method %in% c("glm", "ridge_fallback"))
  expect_true(is.na(fit$ridge_lambda) || is.finite(fit$ridge_lambda))
  p <- predict(fit, d$counts[, 1:6, drop = FALSE], type = "response")
  expect_true(all(is.finite(p)))
})
