test_that("predict works on detectPanel_result objects", {

  counts <- matrix(
    c(
      10,20,30,40,
      20,30,40,50,
      30,40,50,60
    ),
    nrow = 3,
    dimnames = list(
      c("miR-001", "miR-002", "miR-003"),
      c("S1", "S2", "S3", "S4")
    )
  )

  x <- data.frame(
    x1 = c(0, 1, 0, 1),
    x2 = c(1, 0, 1, 0),
    x3 = c(0, 0, 1, 1)
  )
  y <- c(0, 1, 0, 1)
  model <- glm(y ~ x1 + x2 + x3, data = cbind(x, y), family = binomial())

  fake_fit <- list(
    features = rownames(counts),
    transform_method = "log2_cpm",
    prior_count = 1,
    center = c("miR-001"=1,"miR-002"=1,"miR-003"=1),
    scale = c("miR-001"=1,"miR-002"=1,"miR-003"=1),
    model = model,
    predictor_names = c("x1", "x2", "x3"),
    threshold = 0.5,
    positive = "Case",
    negative = "Control"
  )
  class(fake_fit) <- "detectPanel_fit"

  object <- list(final_model = fake_fit)
  class(object) <- "detectPanel_result"

  pred <- predict(object, counts, type = "response")

  expect_equal(length(pred), ncol(counts))
  expect_true(all(is.finite(pred)))
})
