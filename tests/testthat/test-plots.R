test_that("plotting helpers return ggplot objects", {
  d <- make_synthetic_panel_data()
  expr <- transform_assay(d$counts, "log2_cpm")
  p1 <- plot_marker_expression(expr, d$metadata, c("m1", "m2"), "group")
  p2 <- plot_marker_roc(expr, d$y, c("m1", "m2"), positive = "Case")
  de <- data.frame(feature = rownames(expr),
                   log2FoldChange = rnorm(nrow(expr)),
                   padj = runif(nrow(expr)))
  p3 <- plot_volcano(de)
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
})

test_that("confusion metrics use a fixed threshold", {
  m <- confusion_metrics(c(0, 0, 1, 1), c(0.1, 0.7, 0.6, 0.9), 0.5)
  expect_equal(m$TP, 2)
  expect_equal(m$FP, 1)
  expect_equal(m$TN, 1)
})
