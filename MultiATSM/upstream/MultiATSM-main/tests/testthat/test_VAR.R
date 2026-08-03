library(testthat)
library(MultiATSM)

# Load data
LoadData("CM_2024")
RiskFactors <- rbind(GlobalMacro, DomMacro[1:4, ])

# 1) Test unconstrained VAR
test_that("VAR returns correct output for unconstrained case", {
  K <- nrow(RiskFactors)
  res <- VAR(RiskFactors, VARtype = "unconstrained")
  expect_type(res, "list")
  expect_true(all(c("K0Z", "K1Z", "SSZ") %in% names(res)))
  expect_equal(dim(res$K0Z), c(K, 1))
  expect_equal(dim(res$K1Z), c(K, K))
  expect_equal(dim(res$SSZ), c(K, K))
  expect_false(anyNA(res)) # tests whether all entries are numeric
})


# 2) Test constrained VAR
test_that("VAR returns correct output for constrained case", {
  K <- nrow(RiskFactors)
  Bcon_Mat <- matrix(0, nrow = K, ncol = K + 1)
  Bcon_Mat[, 1:3] <- NaN
  res <- VAR(RiskFactors, VARtype = "constrained", Bcon_Mat)
  expect_type(res, "list")
  expect_true(all(c("K0Z", "K1Z", "SSZ") %in% names(res)))
  expect_equal(dim(res$K0Z), c(K, 1))
  expect_equal(dim(res$K1Z), c(K, K))
  expect_equal(dim(res$SSZ), c(K, K))
  expect_false(anyNA(res)) # tests whether all entries are numeric
})

# 3) Test error for invalid VARtype
test_that("VAR throws error for invalid VARtype", {
  expect_error(VAR(RiskFactors, VARtype = "invalid"))
})
