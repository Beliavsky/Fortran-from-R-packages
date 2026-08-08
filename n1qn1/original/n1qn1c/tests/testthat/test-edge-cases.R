fr3 <- function(x) {
  p <- 100
  f <- 1.0
  for (i in 2:3) f <- f + p * (x[i] - x[i-1]^2)^2 + (1.0 - x[i])^2
  f
}
grr3 <- function(x) {
  p <- 100
  g <- double(3)
  g[1] <- -4.0 * p * (x[2] - x[1]^2) * x[1]
  g[2] <- 2.0*p*(x[2]-x[1]^2) - 4.0*p*(x[3]-x[2]^2)*x[2] - 2.0*(1.0-x[2])
  g[3] <- 2.0 * p * (x[3] - x[2]^2) - 2.0 * (1.0 - x[3])
  g
}
x0 <- c(1.02, 1.02, 1.02)

test_that("restart=TRUE sets mode 3 and converges", {
  r1 <- n1qn1(fr3, grr3, x0)
  # restart with previously computed hessian; should still converge
  r2 <- n1qn1(fr3, grr3, x0, zm = r1$c.hess, restart = TRUE)
  expect_equal(r2$par, c(1, 1, 1), tolerance = 1e-3)
  expect_equal(r2$value, 1, tolerance = 1e-6)
})

test_that("assign=TRUE stores c.hess in environment", {
  e <- new.env(parent = emptyenv())
  ret <- n1qn1(fr3, grr3, x0, environment = e, assign = TRUE)
  expect_true(!is.null(e$c.hess))
  expect_identical(e$c.hess, ret$c.hess)
})

test_that("single-variable optimization converges", {
  f1 <- function(x) x^2
  g1 <- function(x) 2 * x
  ret <- n1qn1(f1, g1, c(3.0))
  expect_equal(ret$par, c(0), tolerance = 1e-3)
  expect_equal(ret$value, 0, tolerance = 1e-6)
})

test_that("repeated calls do not crash or corrupt state", {
  results <- lapply(seq_len(5), function(i) n1qn1(fr3, grr3, x0))
  for (r in results) {
    expect_equal(r$par, c(1, 1, 1), tolerance = 1e-3)
  }
})

test_that("error on empty vars (n=0)", {
  expect_error(n1qn1(fr3, grr3, numeric(0)), "n must be positive")
})

test_that("error on wrong-length zm", {
  expect_error(
    n1qn1(fr3, grr3, x0, zm = double(5)),
    "Compressed Hessian not the right length"
  )
})
