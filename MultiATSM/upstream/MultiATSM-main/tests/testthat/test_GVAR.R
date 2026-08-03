library(testthat)
library(MultiATSM)

# Load inputs
data(GVARFactors)

GVARinputs <- list(
  Economies = c("China", "Brazil", "Mexico", "Uruguay"),
  GVARFactors = GVARFactors, VARXtype = "unconstrained"
)

GVARinputs$Wgvar <- matrix(c(
  0, 0.83, 0.86, 0.38,
  0.65, 0, 0.13, 0.55,
  0.32, 0.12, 0, 0.07,
  0.03, 0.05, 0.01, 0
), nrow = 4, ncol = 4)
N <- 3


# 1) Test output structure
test_that("GVAR returns correct output structure", {
  res <- GVAR(GVARinputs, N)
  expect_type(res, "list")
  expect_true(all(c("VARX", "Gy.0", "F0", "F1", "Sigma_y") %in% names(res)))
  # Test for the absence of  NA, NaN, or Inf
  expect_false(any(!is.finite(unlist(res))))
})


# 2) Test error for inconsistent N
test_that("GVAR issues an error for inconsistent N", {
  N <- 10
  expect_error(GVAR(GVARinputs, N, CheckInputs = TRUE))
})

# 3) Test for restricted VARXtype elements
test_that("GVAR checks for a constrained model version (Spanned Factors)", {
  G <- length(GVARFactors$Global)
  GVARinputs$VARXtype <- "constrained: Spanned Factors"
  res <- GVAR(GVARinputs, N, CheckInputs = TRUE)
  expect_true(any(res$F1[-(1:G), ] == 0))
})


# 4) Test for restricted VARXtype elements
test_that("GVAR checks for a constrained model version (Inflation)", {
  GVARinputs$VARXtype <- "constrained: Inflation"
  res <- GVAR(GVARinputs, N, CheckInputs = TRUE)
  infl_rows <- grep("Inflation", rownames(res$F1))
  expect_true(length(infl_rows) > 0)
})

# 5) Test error for missing required GVARinputs elements
test_that("GVAR throws error for invalid (incomplete) VARXtype element", {
  GVARinputs$VARXtype <- "constrained"
  expect_error(GVAR(GVARinputs, N, CheckInputs = TRUE))
})
