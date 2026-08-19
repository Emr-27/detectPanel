test_that("validation aligns samples and encodes outcome", {
  d <- make_synthetic_panel_data()
  shuffled <- d$metadata[sample(rownames(d$metadata)), , drop = FALSE]
  x <- validate_assay_data(d$counts, shuffled, "group", positive = "Case")
  expect_identical(colnames(x$assay), rownames(x$metadata))
  expect_equal(sum(x$outcome), sum(d$y == "Case"))
})

test_that("feature exclusion produces an audit", {
  d <- make_synthetic_panel_data()
  rownames(d$counts)[9:10] <- c("CTRL_a", "HK_b")
  x <- exclude_features(d$counts, regex = "^(CTRL_|HK_)")
  expect_equal(nrow(x$assay), 8)
  expect_setequal(x$excluded$feature, c("CTRL_a", "HK_b"))
})

test_that("log2 CPM is finite", {
  d <- make_synthetic_panel_data()
  x <- transform_assay(d$counts, "log2_cpm")
  expect_equal(dim(x), dim(d$counts))
  expect_true(all(is.finite(x)))
})

test_that("PCA QC flags without deleting samples", {
  d <- make_synthetic_panel_data()
  expr <- transform_assay(d$counts, "log2_cpm")
  qc <- flag_pca_outliers(expr, group = d$y, n_components = 3)
  expect_equal(nrow(qc), ncol(expr))
  expect_true(all(c("sample", "robust_pca_distance", "flagged") %in% names(qc)))
  expect_false(any(qc$flagged))
})
