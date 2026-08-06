# Parametric bootstrap standard errors (added in 1.4-0).

test_that("rsdc_bootstrap returns coherent SE / CI for a tvtp fit", {
  skip_on_cran()
  set.seed(1)
  X <- cbind(1, as.numeric(scale(seq_len(300))))
  y <- scale(matrix(rnorm(300 * 2), 300, 2))
  fit <- rsdc_estimate("tvtp", residuals = y, N = 2, X = X,
                       control = list(itermax = 25, NP = 40, seed = 1))

  bs <- rsdc_bootstrap(fit, B = 20, seed = 1)
  expect_equal(length(bs$se), fit$npar)
  expect_true(all(bs$se >= 0))
  expect_equal(dim(bs$vcov), c(fit$npar, fit$npar))
  expect_equal(nrow(bs$ci), fit$npar)
  expect_true(all(bs$ci[, 1] <= bs$ci[, 2]))
  expect_true(bs$B >= 2L)
  # rho replicates must stay inside (-1, 1) — the bound the Wald interval can violate
  rho_cols <- grep("^rho", colnames(bs$replicates))
  expect_true(all(abs(bs$replicates[, rho_cols]) < 1))
})

test_that("rsdc_bootstrap works for noX (fixed transition matrix)", {
  skip_on_cran()
  set.seed(2)
  y <- scale(matrix(rnorm(250 * 2), 250, 2))
  fit <- rsdc_estimate("noX", residuals = y, N = 2,
                       control = list(itermax = 25, NP = 40, seed = 1))
  bs <- rsdc_bootstrap(fit, B = 15, seed = 2)
  expect_equal(length(bs$se), fit$npar)
  expect_true(all(is.finite(bs$se)))
})

test_that("bootstrap replicates are distinct (RNG regression, 1.6-0 fix)", {
  set.seed(11)
  y <- scale(matrix(rnorm(200 * 2), 200, 2))
  fit <- suppressWarnings(rsdc_estimate("const", residuals = y))
  bs <- suppressWarnings(rsdc_bootstrap(fit, B = 6, seed = 3))
  # before the fix, replicates 2..B were identical: every warm-started refit
  # called set.seed() and so reset the global RNG stream between simulations
  expect_gt(nrow(unique(round(bs$replicates, 10))), 2L)
})
