test_that("detectability and marker AUC identify signal", {
  d <- make_synthetic_panel_data()
  y <- as.integer(d$y == "Case")
  detect <- score_detectability(d$counts, group = y, detection_threshold = 5)
  expr <- transform_assay(d$counts, "log2_cpm")
  auc <- evaluate_markers(expr, y)
  selected <- preselect_markers(
    detect, auc,
    min_mean = 5, min_median = 2,
    min_detection = 0.2, min_group_detection = 0.1,
    min_auc = 0.55, fallback_n = 6
  )
  expect_true(all(c("m1", "m2", "m3") %in% auc$feature[1:5]))
  expect_gte(nrow(selected), 3)
  expect_true(all(is.finite(selected$candidate_score)))
})
