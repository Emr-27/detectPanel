test_that("nested validation returns outer predictions", {
  d <- make_synthetic_panel_data(n = 30, p = 8, seed = 4)
  fit <- nested_validate_panels(
    d$counts, d$metadata, outcome = "group", positive = "Case",
    panel_size = 2, candidate_n = 5,
    detection_threshold = 5,
    min_mean = 5, min_median = 2,
    min_detection = 0.2, min_group_detection = 0.1,
    min_auc = 0.50,
    outer_v = 3, outer_repeats = 1,
    inner_v = 3, inner_repeats = 1,
    seed = 7
  )
  expect_s3_class(fit, "detectPanel_nested")
  expect_gt(nrow(fit$predictions), 0)
  expect_gt(nrow(fit$panel_frequency), 0)
  expect_true(is.finite(fit$outer_summary$pooled_sample_level_AUC))
})
