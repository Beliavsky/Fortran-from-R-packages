if (requireNamespace("testthat", quietly = TRUE)) {
  library(testthat)
  library(PortfolioTesteR)
  test_check("PortfolioTesteR")
} else {
  message("Skipping tests: 'testthat' not installed")
}
