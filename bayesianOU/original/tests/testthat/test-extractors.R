test_that("build_accounting_block returns correct structure", {
  TMG_raw <- rnorm(100)
  zTMG_exo <- (TMG_raw - mean(TMG_raw[1:70])) / sd(TMG_raw[1:70])
  zTMG_use <- zTMG_exo
  mu_tmg <- mean(TMG_raw[1:70])
  sd_tmg <- sd(TMG_raw[1:70])
  
  result <- build_accounting_block(
    TMG_raw, zTMG_exo, zTMG_use,
    mu_tmg, sd_tmg,
    hard = TRUE,
    sigma_delta = 0.002
  )
  
  expect_type(result, "list")
  expect_named(result, c("tmg_byK", "tmg_exo", "wedge_delta",
                         "sigma_delta_prior", "backtransform_max_err", "note"))
  expect_length(result$tmg_byK, 100)
  expect_true(all(result$wedge_delta == 0))
  # Back-transform of the exogenous TMG must reproduce the raw series.
  expect_lt(result$backtransform_max_err, 1e-8)
})

test_that("build_accounting_block soft constraint has non-zero wedge", {
  TMG_raw <- rnorm(100)
  zTMG_exo <- (TMG_raw - mean(TMG_raw[1:70])) / sd(TMG_raw[1:70])
  zTMG_use <- zTMG_exo + rnorm(100, 0, 0.1)
  mu_tmg <- mean(TMG_raw[1:70])
  sd_tmg <- sd(TMG_raw[1:70])
  
  result <- build_accounting_block(
    TMG_raw, zTMG_exo, zTMG_use,
    mu_tmg, sd_tmg,
    hard = FALSE,
    sigma_delta = 0.002
  )
  
  expect_false(all(result$wedge_delta == 0))
})

test_that("drift_decomposition_grid returns correct structure", {
  summ <- list(
    kappa_s = c(0.1, 0.2, 0.3),
    a3_s = c(-0.05, -0.1, -0.15)
  )
  
  result <- drift_decomposition_grid(NULL, summ)
  
  expect_type(result, "list")
  expect_named(result, c("z", "drift"))
  expect_length(result$z, 101)
  expect_equal(ncol(result$drift), 3)
  expect_equal(nrow(result$drift), 101)
})

test_that("drift_decomposition_grid custom grid works", {
  summ <- list(
    kappa_s = c(0.1, 0.2),
    a3_s = c(-0.05, -0.1)
  )
  
  z_grid <- seq(-1, 1, length.out = 21)
  result <- drift_decomposition_grid(NULL, summ, z_grid = z_grid)
  
  expect_length(result$z, 21)
  expect_equal(nrow(result$drift), 21)
})

test_that("drift has correct sign at equilibrium", {
  summ <- list(
    kappa_s = c(0.2),
    a3_s = c(-0.1)
  )
  
  result <- drift_decomposition_grid(NULL, summ, z_grid = seq(-2, 2, by = 0.5))
  
  negative_z_drift <- result$drift[result$z < 0, 1]
  positive_z_drift <- result$drift[result$z > 0, 1]
  
  expect_true(all(negative_z_drift > 0))
  expect_true(all(positive_z_drift < 0))
})
