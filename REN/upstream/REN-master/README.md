‘REN’: Regularization Ensemble for Portfolio Optimization
================

<!-- badges: start -->

[![R-CMD-check](https://github.com/bonsook/REN/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bonsook/REN/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

## Overview

The ‘REN’ package provides a set of tools for performing portfolio
optimization using various regularization and ensemble learning methods.
The package is designed to produce stable out-of-sample return
predictions, particularly in the presence of strong correlations between
assets. The core functions enable users to prepare data, set up parallel
processing, and perform in-depth portfolio analysis.

## Installation

To install the ‘REN’ package, you can use the following command:

``` r
# Install the development version from GitHub
devtools::install_github("bonsook/REN")
```

## Main Functions

### 1. `setup_parallel()`

This function sets up parallel processing to speed up the computation of
portfolio optimization tasks.

- **Parameters:**
  - `num_cores`: The number of cores to use for parallel processing.
    Default is 7.
- **Example:**

``` r
# Set up parallel processing with the default number of cores
cl <- setup_parallel()

# Set up parallel processing with 4 cores
cl <- setup_parallel(num_cores = 4)

# Stop the cluster after completing the analysis
stopCluster(cl) 
```

### 2. `prepare_data()`

This function prepares the input data for portfolio optimization by
structuring it into the required format and calculating necessary
metrics.

- **Parameters:**
  - `data`: A data frame containing the asset returns and other relevant
    metrics.
- **Example:**

``` r
# Prepare the data for analysis
ff25 <- read.csv("data/FF25.csv")

# Define the date column index, start date, and end date
date_column_index <- 1  # Update this based on your data
start_date <- "19990101"  # Adjust as needed
end_date <- "20231231"  # Adjust as needed

# Prepare the data for analysis
data_prep <- prepare_data(ff25, date_column_index, start_date, end_date)
x <- data_prep$x
mon <- data_prep$mon
count <- data_prep$count
Date <- data_prep$Date
```

### 3. `perform_analysis()`

This function performs portfolio analysis using various methods such as
Mean-Variance (MV), James-Stein (JM), LASSO, Ridge Regression, and Equal
Weighting (EW). It calculates portfolio weights, turnover, returns,
Sharpe ratios, volatility, and maximum drawdown for each method.

- **Parameters:**
  - `x`: A numeric matrix where each column represents asset returns and
    rows represent time periods.
  - `mon`: A numeric vector representing the number of months since the
    start date for each time period.
  - `count`: A numeric vector indicating the number of entries per
    month.
  - `Date`: A vector of Date objects representing the dates of the time
    periods.
  - `num_cores`: The number of cores to use for parallel processing.
    Default is 7.
- **Returns:**
  - A list containing components such as cumulative return plots,
    cumulative turnover plots, and performance metrics (Sharpe ratio,
    volatility, max drawdown).
- **Example:**

``` r
# Perform the portfolio analysis
result <- perform_analysis(x, mon, count, Date, num_cores)

# Accessing the results
cumulative_return_plot <- result$cumulative_return_plot
turnover_mean <- result$turnover_mean
sharpe_ratio <- result$sharpe_ratio
volatility <- result$volatility
max_drawdown <- result$max_drawdown

# Display the cumulative return plot
print(cumulative_return_plot)
```

## Example Workflow

Here’s an example workflow using the ‘REN’ package:

``` r
# Step 1: Set up parallel processing
cl <- setup_parallel(num_cores = 4)

# Step 2: Prepare the data
data_prep <- prepare_data(ff25, date_column_index, start_date, end_date)
x <- data_prep$x
mon <- data_prep$mon
count <- data_prep$count
Date <- data_prep$Date

# Step 3: Perform portfolio analysis
result <- perform_analysis(x, mon, count, Date, num_cores)

# Step 4: Plot and interpret the results
print(results$cumulative_return_plot)
print(results$turnover_mean)
print(results$sharpe_ratio)

# Remember to stop the cluster after completing the analysis
stopCluster(cl)
```

## Contributing

Contributions to the ‘REN’ package are welcome. Please feel free to
submit a pull request or report any issues you encounter.

## License

This package is licensed under the MIT License.

------------------------------------------------------------------------
