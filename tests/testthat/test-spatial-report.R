test_that("the spatial report extracts DHARMa's Moran statistic", {
  helpers <- new.env(parent = globalenv())
  sys.source(testthat::test_path("..", "..", "inst", "analysis",
                                 "spatial_analysis_helpers.R"), envir = helpers)
  test <- list(statistic = c(observed = 0.024, expected = -0.001, sd = 0.02),
               p.value = 0.003)
  expect_equal(helpers$moran_test_values(test),
               c(observed = 0.024, expected = -0.001, p_value = 0.003))
  expect_error(helpers$moran_test_values(list(estimate = test$statistic,
                                             p.value = test$p.value)),
               "missing observed/expected")
  expect_error(helpers$moran_test_values(list(statistic = test$statistic,
                                             p.value = numeric())),
               "three finite scalar")
})
