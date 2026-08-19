test_that("stratified splits are reproducible and contain both classes", {
  y <- rep(0:1, each = 12)
  a <- make_repeated_stratified_splits(y, v = 3, repeats = 2, seed = 9)
  b <- make_repeated_stratified_splits(y, v = 3, repeats = 2, seed = 9)
  expect_identical(a, b)
  expect_equal(length(a), 6)
  expect_true(all(vapply(a, function(s) "repeat_id" %in% names(s), logical(1))))
  for (s in a) {
    expect_equal(sort(unique(y[s$assessment])), 0:1)
  }
})

test_that("panel search uses a shared split object", {
  d <- make_synthetic_panel_data(n = 30, p = 8)
  y <- as.integer(d$y == "Case")
  expr <- transform_assay(d$counts, "log2_cpm")
  splits <- make_repeated_stratified_splits(y, v = 3, repeats = 1, seed = 3)
  search <- search_panels(
    expr, y, candidates = rownames(expr)[1:6], panel_size = 3,
    splits = splits, near_best_tolerance = 0.01
  )
  expect_s3_class(search, "detectPanel_search")
  expect_length(search$best_panel, 3)
  expect_identical(search$splits, splits)
  expect_true(all(search$results$mean_CV_AUC >= 0 &
                  search$results$mean_CV_AUC <= 1))
})
