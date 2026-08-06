# The `cores` option of rsdc_bootstrap() and rsdc_corr_bands() (added in
# 1.6-0): parallel execution must reproduce the sequential result exactly.

test_that("bootstrap and corr bands are core-count invariant", {
  # Genuine two-regime data so the fit is interior (bands need valid draws)
  K <- 2
  Sig <- array(c(1, 0.1, 0.1, 1,  1, 0.7, 0.7, 1), c(K, K, 2))
  sim <- rsdc_simulate(n = 250, X = matrix(1, 250, 1),
                       beta = matrix(qlogis(0.9), 2, 1),
                       mu = matrix(0, 2, K), sigma = Sig, N = 2, seed = 7)
  fit <- suppressWarnings(rsdc_estimate("noX", residuals = sim$observations,
                                        N = 2))

  b1 <- suppressWarnings(rsdc_bootstrap(fit, B = 4, seed = 5, cores = 1))
  b2 <- suppressWarnings(rsdc_bootstrap(fit, B = 4, seed = 5, cores = 2))
  expect_equal(b1$replicates, b2$replicates)

  c1 <- suppressWarnings(rsdc_corr_bands(fit, B = 20, seed = 5, cores = 1))
  c2 <- suppressWarnings(rsdc_corr_bands(fit, B = 20, seed = 5, cores = 2))
  expect_equal(c1, c2)
})

test_that(".rsdc_lapply validates cores", {
  expect_error(RSDC:::.rsdc_lapply(1:3, identity, cores = 0), "cores")
  expect_error(RSDC:::.rsdc_lapply(1:3, identity, cores = c(1, 2)), "cores")
  expect_identical(RSDC:::.rsdc_lapply(1:3, function(i) i^2, cores = 1),
                   list(1, 4, 9))
})
