# test-perform_analysis.R

library(testthat)
library(ggplot2)
library(doParallel)

test_that("perform_analysis runs without errors and produces expected output", {
  # Set up the test data
  set.seed(123)
  x <- matrix(runif(700), ncol = 10)  # 10 columns (assets), 70 rows (observations)
  mon <- rep(1:10, each = 7)          # Example month identifiers, 7 observations per month
  count <- rep(7, 10)                 # Example count per month (7 entries per month)
  Date <- as.Date('2020-01-01') + 0:69 # Example date sequence (70 days)

  # Set up parallel backend with 2 cores
  num_cores <- 2
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)

  # Run the analysis
  result <- perform_analysis(x, mon, count, Date, num_cores = num_cores)

  # Test that the output is a list
  expect_type(result, "list")

  # Test that the result contains the expected components
  expect_true(all(c("cumulative_return_plot", "cumulative_turnover_plot",
                    "turnover_mean", "sharpe_ratio", "volatility",
                    "max_drawdown", "vw_to_mean", "vw_sharpe_ratio",
                    "vw_volatility", "vw_max_drawdown") %in% names(result)))

  # Test that the plots are ggplot objects
  expect_s3_class(result$cumulative_return_plot, "ggplot")
  expect_s3_class(result$cumulative_turnover_plot, "ggplot")

  # Check that the numeric components are not NULL
  expect_false(is.null(result$turnover_mean))
  expect_false(is.null(result$sharpe_ratio))
  expect_false(is.null(result$volatility))
  expect_false(is.null(result$max_drawdown))
  expect_false(is.null(result$vw_to_mean))
  expect_false(is.null(result$vw_sharpe_ratio))
  expect_false(is.null(result$vw_volatility))
  expect_false(is.null(result$vw_max_drawdown))

  # Stop the parallel cluster after the test
  stopCluster(cl)
})

