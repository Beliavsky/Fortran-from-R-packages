test_that("VARMA simulation accepts Randompack RNG objects", {
  model <- varmapack_model(A = matrix(0.4, 1, 1), B = matrix(0.2, 1, 1),
                            Sig = matrix(1, 1, 1))
  rng1 <- randompack::randompack_rng()
  rng2 <- randompack::randompack_rng()
  rng1$seed(123)
  rng2$seed(123)
  X1 <- model$sim(20, nrep = 3, rng = rng1)
  X2 <- model$sim(20, nrep = 3, rng = rng2)
  tail1 <- rng1$unif(5)
  tail2 <- rng2$unif(5)
  expect_identical(X1, X2)
  expect_identical(tail1, tail2)
  expect_equal(dim(X1), c(1L, 20L, 3L))
})

test_that("white-noise simulation returns its shocks", {
  model <- varmapack_model(Sig = diag(2))
  out <- model$sim(10, nrep = 2, return_shocks = TRUE)
  expect_named(out, c("X", "E"))
  expect_equal(dim(out$X), c(2L, 10L, 2L))
  expect_equal(dim(out$E), c(2L, 10L, 2L))
  expect_equal(out$X, out$E)
})

test_that("time-dependent VARMA means are applied and extended", {
  mu <- matrix(c(10, 20, 11, 21, 12, 22), 2, 3)
  model <- varmapack_model(Sig = diag(2), mu = mu)
  out <- model$sim(5, nrep = 2, return_shocks = TRUE)
  expected <- cbind(mu, mu[, 3], mu[, 3])
  expect_equal(out$X[, , 1] - out$E[, , 1], expected)
  expect_equal(out$X[, , 2] - out$E[, , 2], expected)
})

test_that("simulation preserves shared and replicate-specific startup paths", {
  model <- varmapack_model(A = matrix(0.4, 1, 1), B = matrix(0.2, 1, 1),
                           Sig = matrix(1, 1, 1))
  shared <- matrix(c(1, 2), 1, 2)
  X <- model$sim(5, nrep = 2, X0 = shared)
  expect_equal(drop(X[, 1:2, 1]), drop(shared))
  expect_equal(drop(X[, 1:2, 2]), drop(shared))
  separate <- array(c(1, 2, 3, 4), c(1L, 2L, 2L))
  X <- model$sim(5, nrep = 2, X0 = separate)
  expect_equal(X[, 1:2, , drop = FALSE], separate)
})

test_that("singular innovation covariance is supported", {
  Sig <- matrix(1, 2, 2)
  model <- varmapack_model(Sig = Sig)
  out <- model$sim(6, nrep = 2, return_shocks = TRUE)
  expect_equal(out$X, out$E)
  expect_equal(out$X[1, , ], out$X[2, , ])
  Theta <- model$irf(0)[, , 1]
  expect_equal(Theta %*% t(Theta), Sig, tolerance = 1e-12)
})
