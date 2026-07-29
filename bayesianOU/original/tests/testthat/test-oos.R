# Regression tests for evaluate_oos index handling.
# Bug: for horizons longer than the test window, (T_train+1):(Tn-hh+1) produced
# a descending sequence and indexed out of range. It must now return NA / n_obs 0.

make_summ <- function(S) {
  list(
    theta_s = rep(0, S),
    kappa_s = rep(0.3, S),
    a3_s    = rep(-0.05, S),
    beta0_s = rep(0.2, S),
    beta1   = 0.3,
    gamma   = 0.1
  )
}

test_that("evaluate_oos returns NA / n_obs 0 when horizon exceeds the test window", {
  set.seed(7)
  Tn <- 10L; S <- 2L; T_train <- 7L
  Yz <- matrix(stats::rnorm(Tn * S), Tn, S)
  Xz <- matrix(stats::rnorm(Tn * S), Tn, S)
  COM <- matrix(abs(stats::rnorm(Tn * S)) + 0.5, Tn, S)
  K   <- matrix(1, Tn, S)

  res <- evaluate_oos(make_summ(S), Yz, Xz, stats::rnorm(Tn), T_train,
                      COM_ts = COM, K_ts = K, com_in_mean = TRUE,
                      horizons = c(1L, 8L))

  # h = 8: last origin = 10 - 8 + 1 = 3 < T_train + 1 = 8 -> nothing to evaluate
  expect_true(is.na(res$h8$RMSE))
  expect_equal(res$h8$n_obs, 0L)

  # h = 1: should evaluate (origins 8, 9, 10)
  expect_true(is.finite(res$h1$RMSE))
  expect_gt(res$h1$n_obs, 0L)
})

test_that("evaluate_oos does not error on plausible small inputs", {
  set.seed(8)
  Tn <- 20L; S <- 3L; T_train <- 14L
  Yz <- matrix(stats::rnorm(Tn * S), Tn, S)
  Xz <- matrix(stats::rnorm(Tn * S), Tn, S)
  COM <- matrix(abs(stats::rnorm(Tn * S)) + 0.5, Tn, S)
  K   <- matrix(stats::runif(Tn * S, 1, 5), Tn, S)

  expect_silent(
    evaluate_oos(make_summ(S), Yz, Xz, stats::rnorm(Tn), T_train,
                 COM_ts = COM, K_ts = K, com_in_mean = FALSE,
                 horizons = c(1L, 4L, 8L))
  )
})
