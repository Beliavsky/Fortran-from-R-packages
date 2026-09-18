test_that("mvn returns an n by d matrix", {
  rng <- randompack_rng("pcg64")
  Sigma <- diag(2)
  x <- rng$mvn(5, Sigma)
  expect_true(is.matrix(x))
  expect_equal(dim(x), c(5, 2))
  expect_true(all(is.finite(x)))
})

test_that("mvn accepts a mean vector", {
  rng <- randompack_rng("pcg64")
  Sigma <- diag(3)
  x <- rng$mvn(4, Sigma, mu = c(1, 2, 3))
  expect_true(is.matrix(x))
  expect_equal(dim(x), c(4, 3))
  expect_true(all(is.finite(x)))
})

test_that("mvn accepts mu = NULL", {
  rng <- randompack_rng("pcg64")
  Sigma <- diag(2)
  x <- rng$mvn(3, Sigma, mu = NULL)
  expect_true(is.matrix(x))
  expect_equal(dim(x), c(3, 2))
  expect_true(all(is.finite(x)))
})

test_that("mvn validates Sigma shape and symmetry", {
  rng <- randompack_rng("pcg64")
  Sigma_ns <- matrix(1, 2, 3)
  expect_error(rng$mvn(2, Sigma_ns))
  Sigma_as <- matrix(c(1, 2, 3, 4), 2, 2)
  expect_error(rng$mvn(2, Sigma_as))
})

test_that("mvn validates mu length", {
  rng <- randompack_rng("pcg64")
  Sigma <- diag(2)
  expect_error(rng$mvn(2, Sigma, mu = c(1, 2, 3)))
})
