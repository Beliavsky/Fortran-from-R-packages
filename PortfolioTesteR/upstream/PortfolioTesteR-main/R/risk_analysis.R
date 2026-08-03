# risk_analysis.R
# Risk analysis and regime detection functions

#library(data.table)
#' Convert Continuous Indicator to Discrete Regimes
#'
#' @description
#' Transforms continuous indicators into discrete regime categories.
#'
#' @param indicator Numeric vector or data frame with indicator values
#' @param breakpoints Numeric vector of breakpoints
#' @param labels Optional character vector of regime names
#' @param use_percentiles Use percentiles instead of fixed breakpoints (default: FALSE)
#'
#' @return Integer vector of regime classifications
#' @export
#' @examples
#' data("sample_prices_weekly")
#' momentum <- calc_momentum(sample_prices_weekly, lookback = 12)
#' selected <- filter_top_n(momentum, 10)
#' # Create VIX-like indicator from volatility
#' vol <- calc_rolling_volatility(sample_prices_weekly, lookback = 20)
#' vix_proxy <- vol$SPY * 100  # Scale to VIX-like values
#' regimes <- create_regime_buckets(vix_proxy, c(15, 25))
create_regime_buckets <- function(indicator, breakpoints, labels = NULL, use_percentiles = FALSE) {
  # Convert continuous indicator into discrete regime buckets
  #
  # Args:
  #   indicator: Numeric vector or single-column data.frame with values to bucket
  #   breakpoints: Numeric vector of breakpoints
  #   labels: Character vector of labels for each bucket (optional)
  #   use_percentiles: If TRUE, breakpoints are percentiles (0-100)
  #                   If FALSE, breakpoints are absolute values (default)

  # Extract values if input is a data structure
  if (is.data.frame(indicator) || is.data.table(indicator)) {
    if (ncol(indicator) > 1) {
      stop("create_regime_buckets: indicator must be a single column")
    }
    values <- indicator[[1]]
  } else {
    values <- indicator
  }

  # Remove NAs for calculations
  valid_values <- values[!is.na(values)]

  if (length(valid_values) == 0) {
    stop("create_regime_buckets: no valid values in indicator")
  }

  # Process breakpoints based on type
  if (use_percentiles) {
    # Treat as percentiles
    if (any(breakpoints <= 0 | breakpoints >= 100)) {
      stop("Percentile breakpoints must be between 0 and 100")
    }
    actual_breaks <- quantile(valid_values, probs = breakpoints / 100)
    actual_breaks <- c(min(valid_values) - 1e-10, actual_breaks, max(valid_values) + 1e-10)
  } else {
    # Treat as absolute values
    actual_breaks <- c(-Inf, breakpoints, Inf)
  }

  # Create labels if not provided
  n_buckets <- length(actual_breaks) - 1
  if (is.null(labels)) {
    if (n_buckets == 2) {
      labels <- c("low", "high")
    } else if (n_buckets == 3) {
      labels <- c("low", "normal", "high")
    } else if (n_buckets == 4) {
      labels <- c("low", "normal", "elevated", "extreme")
    } else {
      labels <- paste0("bucket_", 1:n_buckets)
    }
  }

  # Check label count
  if (length(labels) != n_buckets) {
    stop(sprintf("create_regime_buckets: need %d labels but got %d",
                 n_buckets, length(labels)))
  }

  # Cut the data
  regimes <- cut(values, breaks = actual_breaks, labels = labels, include.lowest = TRUE)

  return(regimes)
}
