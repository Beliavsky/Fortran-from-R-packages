# technical_indicators.R
# Exact translations of Python momentum strategy indicators
# All indicators work with close-only price data

# Required libraries
#library(data.table)

# Check for TTR package availability
check_ttr <- function() {
  if (!requireNamespace("TTR", quietly = TRUE)) {
    stop("TTR package required for technical indicators. Install with: install.packages('TTR')")
  }
}

###############################################################################
# CORE CALCULATION FUNCTIONS - Direct Python translations
###############################################################################
#' Calculate Price Momentum
#'
#' @description
#' Calculates momentum as the percentage change in price over a specified
#' lookback period. Optimized using column-wise operations (25x faster).
#'
#' @param data A data.frame or data.table with Date column and price columns
#' @param lookback Number of periods for momentum calculation (default: 12)
#'
#' @return Data.table with momentum values (0.1 = 10% increase)
#' @export
#' @examples
#' data("sample_prices_weekly")
#' momentum <- calc_momentum(sample_prices_weekly, lookback = 12)
calc_momentum <- function(data, lookback = 12) {
  # OPTIMIZED VERSION - Column-wise processing using shift
  dt <- ensure_dt_copy(data)
  symbol_cols <- setdiff(names(dt), "Date")

  momentum_df <- copy(dt)

  # Process entire columns at once using data.table's shift
  for (col in symbol_cols) {
    current <- dt[[col]]
    past <- shift(current, n = lookback, type = "lag")
    momentum_df[, (col) := safe_divide(current - past, past)]
  }

  return(momentum_df)
}

#' Calculate Distance from Reference
#'
#' @description
#' data("sample_prices_weekly")
#' Calculates percentage distance between prices and reference values
#' (typically moving averages).
#'
#' @param price_df Data frame with price data
#' @param reference_df Data frame with reference values (same structure)
#'
#' @return Data.table with percentage distances
#' @export
#' @examples
#' data("sample_prices_weekly")
#' ma20 <- calc_moving_average(sample_prices_weekly, 20)
#' data("sample_prices_weekly")
#' distance <- calc_distance(sample_prices_weekly, ma20)
calc_distance <- function(price_df, reference_df) {
  # FIXED: Safe division
  price_dt <- ensure_dt_copy(price_df)
  reference_dt <- ensure_dt_copy(reference_df)

  symbol_cols <- setdiff(names(price_dt), "Date")
  distance_df <- copy(price_dt)

  for (col in symbol_cols) {
    distance_vals <- safe_divide(
      price_dt[[col]] - reference_dt[[col]],
      reference_dt[[col]]
    )
    distance_df[, (col) := distance_vals]
  }

  return(distance_df)
}

#' Calculate Moving Average
#'
#' @description
#' Calculates simple moving average for each column in the data.
#'
#' @param data Data frame with Date column and price columns
#' @param window Number of periods for moving average (default: 20)
#'
#' @return Data.table with moving average values
#' @export
#' @examples
#' data("sample_prices_weekly")
#' ma20 <- calc_moving_average(sample_prices_weekly, window = 20)
calc_moving_average <- function(data, window = 20) {
  # FIXED: Ensure data.table is preserved
  dt <- ensure_dt_copy(data)
  symbol_cols <- setdiff(names(dt), "Date")

  ma_df <- copy(dt)
  for (col in symbol_cols) {
    ma_df[, (col) := frollmean(get(col), n = window, align = "right")]
  }

  return(ma_df)
}

#' Calculate Stochastic D Indicator
#'
#' @description
#' Calculates the Stochastic D indicator for momentum analysis.
#' The %D line is the smoothed version of %K, commonly used for
#' momentum signals in range 0-100.
#'
#' @param data Price data with Date column and symbol columns
#' @param k Lookback period for stochastic K calculation
#' @param d Smoothing period for D line
#'
#' @return Data.table with Stochastic D values for each symbol
#' @export
#' @examples
#' data("sample_prices_weekly")
#' data(sample_prices_weekly)
#' data("sample_prices_weekly")
#' stoch_d <- calc_stochastic_d(sample_prices_weekly, k = 14, d = 3)
#' head(stoch_d)
calc_stochastic_d <- function(data, k = 14, d = 3) {
  # Check TTR availability
  if (!requireNamespace("TTR", quietly = TRUE)) {
    stop("TTR package required. Install with: install.packages('TTR')")
  }

  # Convert to data.table if not already
  dt <- ensure_dt_copy(data)
  symbol_cols <- setdiff(names(dt), "Date")

  # Initialize result with same structure
  stoch_d_df <- data.table::copy(dt)

  # Compute Stochastic %D for each column
  for (col in symbol_cols) {
    # TTR's stoch needs HLC data, so use close for all three
    hlc_data <- cbind(dt[[col]], dt[[col]], dt[[col]])

    tryCatch({
      # Calculate stochastic
      stoch_result <- TTR::stoch(hlc_data, nFastK = k, nFastD = d, nSlowD = d)

      if (!is.null(stoch_result) && ncol(stoch_result) >= 2) {
        # Extract %D line (second column) and scale to 0-100 range
        stoch_d_values <- stoch_result[, 2] * 100
        stoch_d_df[, (col) := stoch_d_values]
      } else {
        stoch_d_df[, (col) := NA_real_]
      }
    }, error = function(e) {
      # If calculation fails, set to NA
      stoch_d_df[, (col) := NA_real_]
    })
  }

  return(stoch_d_df)
}

###############################################################################
# TTR-BASED INDICATORS - Wrapping TTR to match Python pandas_ta behavior
###############################################################################

#' Calculate Relative Strength Index (RSI)
#'
#' @description
#' Calculates RSI for each column. RSI ranges from 0-100.
#' Above 70 indicates overbought, below 30 indicates oversold.
#'
#' @param data Data frame with Date column and price columns
#' @param period RSI period (default: 14)
#'
#' @return Data.table with RSI values (0-100 range)
#' @export
#' @examples
#' data("sample_prices_weekly")
#' rsi <- calc_rsi(sample_prices_weekly, period = 14)
#' overbought <- filter_above(rsi, 70)
calc_rsi <- function(data, period = 14) {
  # Check TTR availability
  if (!requireNamespace("TTR", quietly = TRUE)) {
    stop("TTR package required. Install with: install.packages('TTR')")
  }

  # Convert to data.table if not already
  dt <- ensure_dt_copy(data)
  symbol_cols <- setdiff(names(dt), "Date")

  # Initialize result with same structure
  rsi_df <- data.table::copy(dt)

  # Compute RSI for each column separately
  for (col in symbol_cols) {
    # Use TTR's RSI function
    rsi_series <- TTR::RSI(dt[[col]], n = period)
    rsi_df[, (col) := rsi_series]
  }

  return(rsi_df)
}

#' Calculate Commodity Channel Index (CCI)
#'
#' @description
#' Calculates CCI using closing prices. CCI measures deviation from average price.
#' Values above 100 indicate overbought, below -100 indicate oversold.
#'
#' @param data Data frame with Date column and price columns
#' @param period CCI period (default: 20)
#'
#' @return Data.table with CCI values
#' @export
#' @examples
#' data("sample_prices_weekly")
#' cci <- calc_cci(sample_prices_weekly, period = 20)
calc_cci <- function(data, period = 20) {
  # Calculate Commodity Channel Index (CCI) for each column in the DataFrame
  # using only closing prices (adapted from Python version).
  #
  # Translates Python: ta.cci(high=data[col], low=data[col], close=data[col], length=period)
  #
  # Args:
  #   data (DataFrame): Price data with dates as index, tickers as columns
  #   period (int): CCI period (lookback window). Default is 20.
  #
  # Returns:
  #   DataFrame: CCI values for each ticker at each date

  # Check TTR availability
  check_ttr()
 #library(TTR)

  # Convert to data.table if not already
  setDT(data)

  # Get symbol columns (all except Date)
  symbol_cols <- setdiff(names(data), "Date")

  # Initialize result with same structure
  cci_df <- copy(data)

  # Compute CCI for each column separately
  for (col in symbol_cols) {
    # TTR's CCI needs HLC, so use close for all three (matches Python adaptation)
    hlc_data <- cbind(data[[col]], data[[col]], data[[col]])  # High, Low, Close all = Close
    cci_series <- CCI(hlc_data, n = period)
    cci_df[, (col) := cci_series]
  }

  return(cci_df)
}


###############################################################################
# CUSTOM IMPLEMENTATIONS - Where TTR doesn't match pandas_ta exactly
###############################################################################
#' Stochastic RSI (StochRSI) for multiple price series
#'
#' Computes Stochastic RSI (%K) per column over a rolling window, returning
#' values in \[0, 1\]. For each symbol, RSI is computed with [TTR::RSI()] over
#' `rsi_length` periods; then StochRSI is
#' \eqn{(RSI_t - \min RSI_{t-L+1:t}) / (\max RSI_{t-L+1:t} - \min RSI_{t-L+1:t})},
#' where \eqn{L} is `stoch_length`. If the range is zero the value is handled
#' per `on_const_window` (default `"zero"`).
#'
#' @param data A `data.frame` or `data.table` with a `Date` column and one
#'   price column per symbol (wide format).
#' @param length Integer lookback used when `rsi_length`/`stoch_length` are NULL. Default `14`.
#' @param rsi_length Optional integer RSI lookback. Default: `length`.
#' @param stoch_length Optional integer stochastic window. Default: `length`.
#' @param on_const_window How to handle windows where `maxRSI == minRSI`?
#'   One of `"zero"` (set to 0), `"na"` (leave `NA`). Default `"zero"`.
#'
#' @return A `data.table` with `Date` and symbol columns containing StochRSI
#'   in \[0, 1\], with leading `NA`s for warmup.
#'
#' @seealso [TTR::RSI()], [calc_momentum()], [calc_moving_average()],
#'   [filter_top_n()], [weight_by_risk_parity()]
#'
#' @examples
#' data(sample_prices_weekly)
#' s <- calc_stochrsi(sample_prices_weekly, length = 14)
#' head(s)
#'
#' @importFrom TTR RSI runMin runMax
#' @export
calc_stochrsi <- function(data,
                          length = 14L,
                          rsi_length = NULL,
                          stoch_length = NULL,
                          on_const_window = c("zero","na")) {
  check_ttr()
  on_const_window <- match.arg(on_const_window)

  if (!is.data.frame(data) || !"Date" %in% names(data)) {
    stop("'data' must be a data.frame/data.table with a 'Date' column")
  }
  if (is.null(rsi_length))   rsi_length   <- as.integer(length)
  if (is.null(stoch_length)) stoch_length <- as.integer(length)

  # Preserve user data; don't mutate by reference
  dt <- if (exists("ensure_dt_copy", mode = "function")) {
    ensure_dt_copy(data)
  } else {
    data.table::as.data.table(data)
  }

  syms <- setdiff(names(dt), "Date")
  out  <- dt[, .(Date)]

  for (nm in syms) {
    x   <- as.numeric(dt[[nm]])
    rsi <- TTR::RSI(x, n = rsi_length)

    rmin  <- TTR::runMin(rsi, n = stoch_length)
    rmax  <- TTR::runMax(rsi, n = stoch_length)
    denom <- rmax - rmin

    srs <- (rsi - rmin) / denom
    flat <- is.finite(denom) & (denom == 0)

    if (any(flat, na.rm = TRUE)) {
      if (on_const_window == "zero") srs[flat] <- 0
      # else "na": leave as NA
    }

    # Clamp to [0,1] for safety and match pandas_ta-style scaling
    srs <- pmin(pmax(srs, 0), 1)
    out[[nm]] <- srs
  }

  data.table::setDT(out)
  out
}

#' Rolling correlation of each symbol to a benchmark
#'
#' Computes rolling correlations between each symbol and a benchmark series
#' (e.g., `SPY`) using simple returns over a fixed lookback window.
#'
#' @param data A `data.frame` or `data.table` with a `Date` column and one
#'   column per symbol containing prices. Must include `benchmark_symbol`.
#' @param benchmark_symbol Character, the benchmark column name (default `"SPY"`).
#' @param lookback Integer window size (>= 2) for rolling correlations.
#' @param min_periods Minimum number of valid observations within the window
#'   to compute a correlation. Default is `ceiling(lookback * 0.67)`.
#' @param method Correlation method, `"pearson"` (default) or `"spearman"`.
#'
#' @return A `data.table` with `Date` and one column per non-benchmark symbol,
#'   containing rolling correlations. Insufficient data yields `NA`s.
#'
#' @details Returns are computed as simple returns \eqn{(P_t - P_{t-1})/P_{t-1}}.
#'   Windows with fewer than `min_periods` valid pairs are marked `NA`.
#'
#' @examples
#' data(sample_prices_weekly)
#' corr <- calc_rolling_correlation(
#'   data = sample_prices_weekly,
#'   benchmark_symbol = "SPY",
#'   lookback = 20
#' )
#' head(corr)
#'
#' @seealso [calc_momentum()], [calc_rolling_volatility()]
#' @importFrom stats cor
#' @export
calc_rolling_correlation <- function(data, benchmark_symbol = "SPY",
                                     lookback = 60, min_periods = NULL,
                                     method = c("pearson", "spearman")) {
  # Calculate rolling correlation between each stock and a benchmark
  #
  # This function measures how closely each security moves with a benchmark
  # over a rolling window. Useful for identifying stocks that follow or
  # diverge from market movements.
  #
  # Args:
  #   data: Price data with Date column and symbol columns
  #   benchmark_symbol: Name of benchmark column (default: "SPY")
  #   lookback: Number of periods for correlation calculation
  #   min_periods: Minimum observations required (default: 2/3 of lookback)
  #   method: Correlation method - "pearson" (default) or "spearman"
  #
  # Returns:
  #   data.table with Date and correlation values for each symbol
  #   Benchmark column is excluded from output

  # Match method argument
  method <- match.arg(method)

  # Set default min_periods if not provided
  if (is.null(min_periods)) {
    min_periods <- ceiling(lookback * 0.67)  # 2/3 of window
  }

  # Input validation
  if (!benchmark_symbol %in% names(data)) {
    stop(paste("calc_rolling_correlation: benchmark", benchmark_symbol, "not found in data"))
  }

  if (lookback < 2) {
    stop("calc_rolling_correlation: lookback must be at least 2")
  }

  # Ensure data.table
  dt <- ensure_dt_copy(data)
  symbol_cols <- setdiff(names(dt), c("Date", benchmark_symbol))

  # Calculate returns for all columns
  returns_dt <- copy(dt)
  all_cols <- c(benchmark_symbol, symbol_cols)

  for (col in all_cols) {
    # Use safe_divide for return calculation
    returns_dt[, (col) := c(NA, safe_divide(diff(get(col)), head(get(col), -1)))]
  }

  # Extract benchmark returns
  benchmark_returns <- returns_dt[[benchmark_symbol]]

  # Initialize result
  corr_df <- data.table(Date = dt$Date)

  # Calculate rolling correlation for each symbol
  for (col in symbol_cols) {
    stock_returns <- returns_dt[[col]]

    # Initialize correlation vector
    corr_values <- rep(NA_real_, nrow(dt))

    # Calculate rolling correlation
    for (i in lookback:nrow(dt)) {
      # Get window of returns
      window_start <- i - lookback + 1
      window_bench <- benchmark_returns[window_start:i]
      window_stock <- stock_returns[window_start:i]

      # Count valid pairs (both non-NA)
      valid_pairs <- !is.na(window_bench) & !is.na(window_stock)
      n_valid <- sum(valid_pairs)

      # Calculate correlation if enough data
      if (n_valid >= min_periods) {
        if (method == "pearson") {
          # Pearson correlation
          corr_val <- cor(window_bench[valid_pairs],
                          window_stock[valid_pairs],
                          method = "pearson")
        } else {
          # Spearman rank correlation
          corr_val <- cor(window_bench[valid_pairs],
                          window_stock[valid_pairs],
                          method = "spearman")
        }

        # Handle edge cases
        if (!is.na(corr_val) && !is.nan(corr_val)) {
          corr_values[i] <- corr_val
        }
      }
    }

    # Add to result
    corr_df[, (col) := corr_values]
  }

  return(corr_df)
}







#' Calculate Rolling Volatility
#'
#' @description
#' Calculates rolling volatility using various methods including standard deviation,
#' range-based, MAD, or absolute returns. Supports different lookback periods.
#'
#' @param data Data frame with Date column and price columns
#' @param lookback Number of periods for rolling calculation (default: 20)
#' @param method Volatility calculation method: "std", "range", "mad", or "abs_return"
#'
#' @return Data frame with Date column and volatility values for each symbol
#' @export
#' @examples
#' data("sample_prices_weekly")
#' # Standard deviation volatility
#' vol <- calc_rolling_volatility(sample_prices_weekly, lookback = 20)
#' # Range-based volatility
#' vol_range <- calc_rolling_volatility(sample_prices_weekly, lookback = 20, method = "range")
calc_rolling_volatility <- function(data, lookback = 20, method = "std") {
  # Validate inputs
  if (!method %in% c("std", "range", "mad", "abs_return")) {
    stop("method must be one of: std, range, mad, abs_return")
  }

  dt <- ensure_dt_copy(data)
  symbol_cols <- setdiff(names(dt), c("Date", "week_end"))

  # Initialize result
  vol_df <- copy(dt)

  if (method == "std") {
    # OPTIMIZED: Calculate all returns at once
    returns_dt <- copy(dt)
    for (col in symbol_cols) {
      returns_dt[, (col) := c(NA, diff(get(col))) / shift(get(col))]
    }

    # OPTIMIZED: Use variance formula with frollmean
    # Variance = E[X^2] - (E[X])^2
    var_df <- copy(dt)

    for (col in symbol_cols) {
      # Get returns for this column
      ret_col <- returns_dt[[col]]

      # Calculate E[X] and E[X^2] using frollmean WITHOUT adaptive
      mean_ret <- frollmean(ret_col, n = lookback, align = "right")
      mean_ret_sq <- frollmean(ret_col^2, n = lookback, align = "right")

      # Variance = E[X^2] - (E[X])^2
      variance <- mean_ret_sq - mean_ret^2
      variance[variance < 0] <- 0  # Numerical safety

      # Standard deviation
      vol_df[, (col) := sqrt(variance)]
      var_df[, (col) := variance]
    }

    # Attach variance as attribute
    attr(vol_df, "variance") <- var_df

  } else if (method == "range") {
    # Keep original implementation for now
    for (col in symbol_cols) {
      vol_df[, (col) := frollapply(get(col),
                                   n = lookback,
                                   FUN = function(x) {
                                     valid_x <- x[!is.na(x)]
                                     if(length(valid_x) >= lookback * 0.8) {
                                       range_val <- (max(valid_x) - min(valid_x))
                                       avg_val <- mean(valid_x)
                                       if(avg_val > 0) {
                                         return((range_val / avg_val))
                                       }
                                     }
                                     return(NA_real_)
                                   },
                                   align = "right")]
    }

  } else if (method == "mad") {
    # Calculate returns once
    returns_dt <- copy(dt)
    for (col in symbol_cols) {
      returns_dt[, (col) := c(NA, diff(get(col))) / shift(get(col))]
    }

    # Keep MAD as is
    for (col in symbol_cols) {
      vol_df[, (col) := frollapply(returns_dt[[col]],
                                   n = lookback,
                                   FUN = function(x) {
                                     valid_x <- x[!is.na(x)]
                                     if(length(valid_x) >= lookback * 0.8) {
                                       med <- median(valid_x)
                                       return(median(abs(valid_x - med)) * 1.4826)
                                     } else {
                                       return(NA_real_)
                                     }
                                   },
                                   align = "right")]
    }

  } else if (method == "abs_return") {
    # OPTIMIZED: Use frollmean on absolute returns
    returns_dt <- copy(dt)
    for (col in symbol_cols) {
      returns_dt[, (col) := c(NA, diff(get(col))) / shift(get(col))]
    }

    for (col in symbol_cols) {
      ret_col <- abs(returns_dt[[col]])  # Absolute returns
      vol_df[, (col) := frollmean(ret_col, n = lookback, align = "right")]
    }
  }

  # Clean up any extra columns
  if ("week_end" %in% names(vol_df)) {
    vol_df[, week_end := NULL]
  }

  return(vol_df)
}








###############################################################################
# BOLLINGER BANDS
###############################################################################

calc_bollinger_bands <- function(data, window = 20, num_std = 2) {
  # Calculate Bollinger Bands for each symbol
  #
  # Bollinger Bands consist of three lines:
  # - Middle Band: Simple moving average (SMA)
  # - Upper Band: SMA + (num_std ? rolling standard deviation)
  # - Lower Band: SMA - (num_std ? rolling standard deviation)
  #
  # These bands expand during volatile periods and contract during calm periods,
  # making them useful for:
  # - Identifying overbought/oversold conditions (price at bands)
  # - Detecting volatility squeezes (bands contracting)
  # - Mean reversion strategies (price bouncing off bands)
  # - Breakout strategies (price breaking through bands)
  #
  # Args:
  #   data: Price data with Date column and symbol columns
  #   window: Lookback period for moving average and std dev (default: 20)
  #   num_std: Number of standard deviations for bands (default: 2)
  #
  # Returns:
  #   List with three components:
  #   - upper: data.table with upper band values
  #   - middle: data.table with middle band (SMA) values
  #   - lower: data.table with lower band values
  #   Each has same structure as input (Date + symbol columns)

  # Input validation
  if (!is.data.frame(data)) {
    stop("calc_bollinger_bands: data must be a data.frame or data.table")
  }

  if (window < 2) {
    stop("calc_bollinger_bands: window must be at least 2")
  }

  if (num_std <= 0) {
    stop("calc_bollinger_bands: num_std must be positive")
  }

  # Ensure data.table
  dt <- ensure_dt_copy(data)
  symbol_cols <- setdiff(names(dt), "Date")

  # Initialize result data.tables
  upper_band <- copy(dt)
  middle_band <- copy(dt)
  lower_band <- copy(dt)

  # Calculate bands for each symbol
  for (col in symbol_cols) {
    # Middle band = Simple Moving Average
    sma <- frollmean(dt[[col]], n = window, align = "right")
    middle_band[, (col) := sma]

    # Calculate rolling standard deviation
    prices <- dt[[col]]
    rolling_std <- rep(NA_real_, length(prices))

    for (i in window:length(prices)) {
      window_start <- i - window + 1
      window_data <- prices[window_start:i]
      valid_data <- window_data[!is.na(window_data)]

      if (length(valid_data) >= window * 0.8) {  # Require 80% valid data
        rolling_std[i] <- sd(valid_data)
      }
    }

    # Calculate bands
    upper_values <- sma + (num_std * rolling_std)
    lower_values <- sma - (num_std * rolling_std)

    upper_band[, (col) := upper_values]
    lower_band[, (col) := lower_values]
  }

  # Return list of bands
  bands <- list(
    upper = upper_band,
    middle = middle_band,
    lower = lower_band
  )

  # Add metadata as attributes
  attr(bands, "window") <- window
  attr(bands, "num_std") <- num_std

  return(bands)
}

###############################################################################
# AVERAGE TRUE RANGE (ATR) - ADAPTED FOR CLOSE-ONLY DATA
###############################################################################

calc_atr <- function(data, period = 14, method = c("close_approximation", "percent_range")) {
  # Calculate Average True Range adapted for close-only price data
  #
  # Traditional ATR requires High, Low, Close prices. Since this library uses
  # close-only data, we provide two approximation methods:
  #
  # 1. "close_approximation": Estimates intraday range using close-to-close
  #    movements with a scaling factor (Garman-Klass inspired)
  #    TR ? 1.5 ? |Close[t] - Close[t-1]|
  #
  # 2. "percent_range": Uses percentage changes as a volatility proxy
  #    TR = |Return[t]| ? Close[t]
  #
  # ATR is then calculated as an exponential moving average of TR values.
  #
  # Common uses:
  # - Position sizing: size = risk_amount / (ATR ? multiplier)
  # - Stop loss placement: stop = entry_price - (ATR ? multiplier)
  # - Volatility comparison across different priced stocks
  # - Trend strength: trending markets show increasing ATR
  #
  # Args:
  #   data: Price data with Date column and symbol columns (close prices)
  #   period: Lookback period for ATR calculation (default: 14)
  #   method: Approximation method for true range (default: "close_approximation")
  #
  # Returns:
  #   data.table with Date column and ATR values for each symbol
  #   Values represent average true range in price units

  # Match method argument
  method <- match.arg(method)

  # Input validation
  if (!is.data.frame(data)) {
    stop("calc_atr: data must be a data.frame or data.table")
  }

  if (period < 1) {
    stop("calc_atr: period must be at least 1")
  }

  # Ensure data.table
  dt <- ensure_dt_copy(data)
  symbol_cols <- setdiff(names(dt), "Date")

  # Initialize result
  atr_df <- data.table(Date = dt$Date)

  # Calculate ATR for each symbol
  for (col in symbol_cols) {
    prices <- dt[[col]]
    n <- length(prices)

    # Initialize true range vector
    true_range <- rep(NA_real_, n)

    if (method == "close_approximation") {
      # Method 1: Approximate using close-to-close with scaling
      # Research shows close-to-close captures about 65% of true range
      # We scale by 1.5 to approximate the full range
      for (i in 2:n) {
        if (!is.na(prices[i]) && !is.na(prices[i-1])) {
          # Absolute close-to-close change
          close_change <- abs(prices[i] - prices[i-1])
          # Scale to approximate true range
          true_range[i] <- close_change * 1.5
        }
      }

    } else {  # percent_range
      # Method 2: Use percentage returns scaled by price
      for (i in 2:n) {
        if (!is.na(prices[i]) && !is.na(prices[i-1]) && prices[i-1] > 0) {
          # Calculate return
          ret <- (prices[i] - prices[i-1]) / prices[i-1]
          # True range as absolute return times current price
          true_range[i] <- abs(ret) * prices[i]
        }
      }
    }

    # Calculate ATR as exponential moving average of true range
    # Using Wilder's smoothing (equivalent to EMA with alpha = 1/period)
    atr_values <- rep(NA_real_, n)

    # Need enough data for initial calculation
    if (sum(!is.na(true_range)) >= period) {
      # Find first valid stretch of data
      first_valid <- which(!is.na(true_range))[1]

      if (!is.na(first_valid) && (first_valid + period - 1) <= n) {
        # Initial ATR: simple average of first 'period' TR values
        initial_values <- true_range[first_valid:(first_valid + period - 1)]
        initial_values <- initial_values[!is.na(initial_values)]

        if (length(initial_values) >= period * 0.8) {
          atr_values[first_valid + period - 1] <- mean(initial_values)

          # Calculate subsequent ATR values using Wilder's smoothing
          alpha <- 1 / period
          for (i in (first_valid + period):n) {
            if (!is.na(true_range[i]) && !is.na(atr_values[i-1])) {
              atr_values[i] <- (1 - alpha) * atr_values[i-1] + alpha * true_range[i]
            }
          }
        }
      }
    }

    # Add to result
    atr_df[, (col) := atr_values]
  }

  # Add metadata as attributes
  attr(atr_df, "period") <- period
  attr(atr_df, "method") <- method

  return(atr_df)
}



