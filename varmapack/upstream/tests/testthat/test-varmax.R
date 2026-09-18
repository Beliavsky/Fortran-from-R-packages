test_that("VARMAX simulation accepts vector and replicate-specific exogenous inputs", {
  A <- array(c(0.2, -0.1, 0.1, 0.3), c(2L, 2L, 1L))
  B <- array(c(0.1, 0.05, -0.2, 0.2), c(2L, 2L, 1L))
  C <- array(c(0.3, 0.1, -0.2, 0.4), c(2L, 2L, 1L))
  model <- varmapack_model(A = A, B = B, C = C, Sig = diag(2))
  X0 <- array(seq_len(8)/10, c(2L, 2L, 2L))
  z <- array(seq_len(40)/10, c(2L, 10L, 2L))
  rng1 <- randompack::randompack_rng()
  rng2 <- randompack::randompack_rng()
  rng1$seed(99)
  rng2$seed(99)
  out1 <- model$sim(10, nrep = 2, X0 = X0, z = z, rng = rng1,
                    return_shocks = TRUE)
  out2 <- model$sim(10, nrep = 2, X0 = X0, z = z, rng = rng2,
                    return_shocks = TRUE)
  expect_identical(out1, out2)
  expect_equal(dim(out1$X), c(2L, 10L, 2L))
  expect_equal(out1$X[, 1:2, , drop = FALSE], X0)
  for (j in 1:2) {
    for (t in 3:10) {
      expected <- out1$E[, t, j] + A[, , 1] %*% out1$X[, t - 1, j] +
        B[, , 1] %*% out1$E[, t - 1, j] + C[, , 1] %*% z[, t, j]
      expect_equal(out1$X[, t, j], drop(expected), tolerance = 1e-12)
    }
  }
})

test_that("VARMAX simulation requires its minimum inputs", {
  model <- varmapack_model(C = matrix(0.3, 1, 1), Sig = matrix(1, 1, 1))
  expect_error(model$sim(10, X0 = matrix(0, 1, 1)))
  model_with_ar <- varmapack_model(A = matrix(0.2, 1, 1),
                                   C = matrix(0.3, 1, 1), Sig = matrix(1, 1, 1))
  expect_error(model_with_ar$sim(10, z = matrix(0, 1, 10)))
})

test_that("VARMAX simulation accepts the minimum startup length", {
  model <- varmapack_model(
    A = matrix(0.6, 1, 1), B = matrix(0.2, 1, 1),
    C = array(c(0.4, -0.1), c(1L, 1L, 2L)), Sig = matrix(1, 1, 1))
  z <- matrix(c(1, -1, 2, -2, 3, -3), 1, 6)
  rng <- randompack::randompack_rng()
  rng$seed(109)
  out <- model$sim(6, nrep = 2, X0 = matrix(0.5, 1, 1), z = z, rng = rng,
                   return_shocks = TRUE)
  for (j in 1:2) {
    for (t in 2:6) {
      expected <- 0.6*out$X[1, t - 1, j] + out$E[1, t, j] +
        0.2*out$E[1, t - 1, j] + 0.4*z[1, t] - 0.1*z[1, t - 1]
      expect_equal(out$X[1, t, j], expected, tolerance = 1e-12)
    }
  }
})

test_that("VARMAX simulation accepts a zero startup length", {
  model <- varmapack_model(C = matrix(0.5, 1, 1), Sig = matrix(1, 1, 1))
  z <- matrix(1:5, 1, 5)
  rng <- randompack::randompack_rng()
  rng$seed(110)
  out <- model$sim(5, nrep = 2, z = z, rng = rng, return_shocks = TRUE)
  for (j in 1:2)
    expect_equal(out$X[1, , j], out$E[1, , j] + 0.5*z[1, ], tolerance = 1e-12)
})

test_that("nonstationary VARMA simulation accepts supplied startup values", {
  model <- varmapack_model(A = matrix(1.25, 1, 1), B = matrix(0.5, 1, 1),
                           Sig = matrix(1, 1, 1))
  rng <- randompack::randompack_rng()
  rng$seed(111)
  out <- model$sim(6, X0 = matrix(c(2, 3, 4), 1, 3), rng = rng,
                   return_shocks = TRUE)
  expect_equal(out$X[1, 1:3, 1], c(2, 3, 4))
  for (t in 2:6) {
    expected <- 1.25*out$X[1, t - 1, 1] + out$E[1, t, 1] +
      0.5*out$E[1, t - 1, 1]
    expect_equal(out$X[1, t, 1], expected, tolerance = 1e-12)
  }
})
