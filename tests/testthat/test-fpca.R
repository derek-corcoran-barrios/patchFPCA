test_that("uncentred fit matches prcomp", {
  set.seed(11)
  profiles <- matrix(stats::rexp(40), nrow = 8, ncol = 5)
  species <- rep(c("species_a", "species_b"), each = 4)

  fit <- fit_patch_fpca(
    profiles,
    species,
    log_transform = TRUE,
    balance = FALSE,
    center = FALSE,
    scale. = FALSE
  )
  direct <- stats::prcomp(log1p(profiles), center = FALSE, scale. = FALSE)

  expect_equal(abs(fit$loadings), abs(direct$rotation), tolerance = 1e-10)
  expect_equal(
    fit$component_share,
    direct$sdev^2 / sum(direct$sdev^2),
    tolerance = 1e-10
  )
  expect_identical(fit$component_share_type, "uncentred_inertia")
})

test_that("projection applies the fitted preprocessing", {
  profiles <- matrix(seq_len(30), nrow = 6, ncol = 5)
  fit <- fit_patch_fpca(profiles, balance = FALSE)
  projected <- project_patch_fpca(fit, profiles)
  expect_equal(unname(as.matrix(projected)), unname(as.matrix(fit$scores)))
})

test_that("centering comparison returns both estimands", {
  set.seed(1)
  profiles <- matrix(stats::rexp(100), ncol = 5)
  species <- rep(c("a", "b"), each = 10)
  comparison <- compare_centering(profiles, species)
  expect_identical(comparison$uncentred$component_share_type, "uncentred_inertia")
  expect_identical(comparison$centred$component_share_type, "centred_variance")
  expect_equal(nrow(comparison$summary), 1L)
})
