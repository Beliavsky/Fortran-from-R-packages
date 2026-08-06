test_that("rsdc_simulate generates consistent shapes", {
  skip_if_not_installed("mvtnorm")

  n <- 25; K <- 3; N <- 2; p <- 2
  X <- cbind(1, rnorm(n))
  beta <- matrix(0, nrow = N, ncol = (N - 1) * p)
  # Encourage self-persistence
  beta[1, ] <- c(1.5, 0.0)
  beta[2, ] <- c(1.2, 0.0)

  mu <- rbind(rep(0, K), rep(0, K))
  # two regimes: different correlation levels
  rho <- toy_rho_matrix(K)
  Sig <- array(0, dim = c(K, K, N))
  for (i in 1:N) {
    R <- diag(K)
    R[lower.tri(R)] <- rho[i, ]
    R[upper.tri(R)] <- t(R)[upper.tri(R)]
    Sig[,,i] <- R
  }

  sim <- rsdc_simulate(n = n, X = X, beta = beta, mu = mu, sigma = Sig, N = N, seed = 99)
  expect_equal(length(sim$states), n)
  expect_equal(dim(sim$observations), c(n, K))
  expect_equal(dim(sim$transition_matrices), c(N, N, n))
  expect_true(all(is.na(sim$transition_matrices[,,1])))
})
