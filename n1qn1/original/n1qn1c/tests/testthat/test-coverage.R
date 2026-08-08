# Shared fixture
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

test_that("print.functions=TRUE exercises the fwrap print path", {
  result <- n1qn1(fr3, grr3, x0, print.functions = TRUE)
  expect_equal(result$par, c(1, 1, 1), tolerance = 1e-3)
})

test_that("invisible=1 sets imp=0 and suppresses print.functions", {
  result <- n1qn1(fr3, grr3, x0, invisible = 1)
  expect_equal(result$par, c(1, 1, 1), tolerance = 1e-3)
})

test_that("invisible != 1 (non-NULL) enables print.functions", {
  result <- n1qn1(fr3, grr3, x0, invisible = 2)
  expect_equal(result$par, c(1, 1, 1), tolerance = 1e-3)
})

test_that("max_iterations without nsim auto-sets nsim = max_iterations * 10", {
  result <- n1qn1(fr3, grr3, x0, max_iterations = 50)
  expect_equal(result$par, c(1, 1, 1), tolerance = 1e-3)
})

test_that("EvalCompiled path works with compiled XPtr function pointers", {
  skip_if_not_installed("Rcpp")
  # Compile helpers that return XPtr<funcPtr> (funcPtr = NumericVector(*)(SEXP,SEXP))
  # This exercises the TYPEOF(fSEXP)==EXTPTRSXP branch in RcppExpMod.cpp and the
  # entire EvalCompiled class in evaluate.h.
  src <- '
#include <Rcpp.h>
typedef Rcpp::NumericVector (*funcPtr)(SEXP, SEXP);
static Rcpp::NumericVector sumSqFn(SEXP par, SEXP) {
  Rcpp::NumericVector x(par);
  return Rcpp::NumericVector::create(Rcpp::sum(x * x));
}
static Rcpp::NumericVector sumSqGr(SEXP par, SEXP) {
  return 2.0 * Rcpp::NumericVector(par);
}
// [[Rcpp::export]]
SEXP makeFnPtr() { return Rcpp::XPtr<funcPtr>(new funcPtr(&sumSqFn)); }
// [[Rcpp::export]]
SEXP makeGrPtr() { return Rcpp::XPtr<funcPtr>(new funcPtr(&sumSqGr)); }
'
  tmp <- tempfile(fileext = ".cpp")
  writeLines(src, tmp)
  Rcpp::sourceCpp(tmp)

  result <- n1qn1(makeFnPtr(), makeGrPtr(), c(1.0, 2.0, 3.0))
  expect_equal(result$par, c(0, 0, 0), tolerance = 1e-3)
  expect_equal(result$value, 0, tolerance = 1e-6)
})
