test_that("cosine similarity removes arbitrary PCA sign", {
  x <- c(1, 2, 3, 4, 5)
  expect_equal(cosine_similarity(x, -x), 1)
  expect_equal(cosine_similarity(x, -x, align_sign = FALSE), -1)
})

test_that("mapping and quality invariance use scenario names", {
  make_result <- function(v) list(loadings = cbind(v, rev(v)))
  results <- list(
    broad_High = make_result(1:5),
    broad_Medium = make_result(1:5 + 0.01),
    broad_Low = make_result(1:5 + 0.02),
    narrow_High = make_result(-(1:5)),
    narrow_Medium = make_result(-(1:5 + 0.01)),
    narrow_Low = make_result(-(1:5 + 0.02))
  )

  mapping <- compare_mapping_invariance(results)
  quality <- compare_quality_invariance(results, "broad")
  expect_equal(mapping$cosine_similarity, rep(1, 3), tolerance = 1e-12)
  expect_equal(nrow(quality), 3L)
})

