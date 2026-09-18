test_that("model analysis methods return expected AR(1) quantities", {
  model <- varmapack_model(A = matrix(0.5, 1, 1), Sig = matrix(1, 1, 1))
  expect_equal(model$specrad(), 0.5)
  expect_equal(model$ma_specrad(), 0)
  expect_equal(drop(model$acvf(1)), c(4/3, 2/3), tolerance = 1e-12)
  expect_equal(drop(model$psi(2)), c(1, 0.5, 0.25), tolerance = 1e-12)
  expect_equal(drop(model$irf(2)), c(1, 0.5, 0.25), tolerance = 1e-12)
})

test_that("multivariate impulse responses preserve matrix orientation", {
  A <- matrix(c(0.1, 0.3, 0.2, 0.4), 2, 2)
  B <- matrix(c(0.5, 0.7, 0.6, 0.8), 2, 2)
  Sig <- matrix(c(1, 0.25, 0.25, 2), 2, 2)
  model <- varmapack_model(A = A, B = B, Sig = Sig)
  expected <- array(c(
    1, 0, 0, 1,
    0.6, 1.0, 0.8, 1.2,
    0.26, 0.58, 0.32, 0.72
  ), c(2L, 2L, 3L))
  expect_equal(model$A[, , 1], A)
  expect_equal(model$B[, , 1], B)
  expect_equal(model$psi(2), expected, tolerance = 1e-12)
  expect_gt(model$specrad(), 0)
  expect_gt(model$ma_specrad(), 0)
})

test_that("sample autocovariance uses series in rows", {
  X <- rbind(1:8, 2*(1:8))
  C <- varmapack_autocov(X, maxlag = 1)
  expect_equal(dim(C), c(2L, 2L, 2L))
  expect_equal(C[1, 2, 1], 2*C[1, 1, 1])
  Xc <- X - rowMeans(X)
  expected <- Xc[, 2:8] %*% t(Xc[, 1:7])/8
  expect_equal(C[, , 2], expected)
  expect_error(varmapack_autocov(X, maxlag = 8))
})

test_that("sample autocovariance supports both normalizations", {
  X <- matrix(c(1, 3, 5, 7), 1, 4)
  Xc <- X - rowMeans(X)
  ml <- varmapack_autocov(X, maxlag = 2, norm = "ML")
  corrected <- varmapack_autocov(X, maxlag = 2, norm = "C")
  for (lag in 0:2) {
    current <- (lag + 1):4
    previous <- 1:(4 - lag)
    product <- Xc[, current, drop = FALSE] %*%
      t(Xc[, previous, drop = FALSE])
    expect_equal(drop(ml[, , lag + 1]), drop(product)/4)
    expect_equal(drop(corrected[, , lag + 1]), drop(product)/(4 - lag))
  }
})

test_that("covariance sequences convert to correlations", {
  cov <- array(c(4, 3, 3, 9, 2, 6, -3, 1.5), dim = c(2, 2, 2))
  original <- cov
  expected <- array(c(1, 0.5, 0.5, 1, 0.5, 1, -0.5, 1/6),
                    dim = c(2, 2, 2))
  corr <- varmapack_cov2corr(cov)
  expect_equal(corr, expected)
  expect_identical(cov, original)
  expect_equal(varmapack_cov2corr(cov[, , 1]), expected[, , 1])
})

test_that("theoretical autocovariance rejects VARMAX models", {
  model <- varmapack_model(C = matrix(0.2, 1, 1), Sig = matrix(1, 1, 1))
  expect_error(model$acvf(1))
})
