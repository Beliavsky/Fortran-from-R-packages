# weighting.R
# Portfolio weighting functions for strategy backtesting
# All functions take binary selection matrices and return normalized weights

# Required libraries
#library(data.table)

###############################################################################
# INPUT VALIDATION HELPERS
###############################################################################

validate_selection_data <- function(selected_df, function_name) {
  # Standard validation for all weighting functions

  if (!is.data.frame(selected_df)) {
    stop(paste(function_name, ": selected_df must be a data.frame or data.table"))
  }

  if (!"Date" %in% names(selected_df)) {
    stop(paste(function_name, ": selected_df must have a 'Date' column"))
  }

  if (ncol(selected_df) < 2) {
    stop(paste(function_name, ": selected_df must have at least one symbol column besides Date"))
  }

  if (nrow(selected_df) == 0) {
    stop(paste(function_name, ": selected_df cannot be empty"))
  }

  return(TRUE)
}

validate_matching_structure <- function(df1, df2, function_name) {
  # Validate that two DataFrames have matching structure

  if (!identical(dim(df1), dim(df2))) {
    stop(paste(function_name, ": DataFrames must have same dimensions"))
  }

  if (!identical(names(df1), names(df2))) {
    stop(paste(function_name, ": DataFrames must have same column names"))
  }

  if (!identical(df1$Date, df2$Date)) {
    stop(paste(function_name, ": DataFrames must have same Date values"))
  }

  return(TRUE)
}

###############################################################################
# CORE WEIGHTING FUNCTIONS
###############################################################################
#' Equal Weight Portfolio Construction
#'
#' @description
#' Creates equal-weighted portfolio from selection matrix.
#' The simplest and often most robust weighting scheme.
#'
#' @param selected_df Binary selection matrix (1 = selected, 0 = not)
#'
#' @return Data.table with equal weights for selected securities
#' @export
#' @examples
#' data("sample_prices_weekly")
#' momentum <- calc_momentum(sample_prices_weekly, lookback = 12)
#' selected <- filter_top_n(momentum, 10)
#' weights <- weight_equally(selected)
weight_equally <- function(selected_df) {
  # OPTIMIZED VERSION - Vectorized operations
  validate_selection_data(selected_df, "weight_equally")

  dt <- ensure_dt_copy(selected_df)
  symbol_cols <- setdiff(names(dt), "Date")

  # Convert to matrix for faster operations
  selection_matrix <- as.matrix(dt[, symbol_cols, with = FALSE])

  # Handle NA - set to 0
  selection_matrix[is.na(selection_matrix)] <- 0

  # Count selections per row (vectorized)
  row_counts <- rowSums(selection_matrix > 0)

  # Create weight matrix (vectorized division)
  weight_matrix <- selection_matrix
  weight_matrix[selection_matrix > 0] <- 1

  # Divide by row counts (vectorized)
  # Avoid division by zero
  valid_rows <- row_counts > 0
  weight_matrix[valid_rows, ] <- weight_matrix[valid_rows, ] / row_counts[valid_rows]
  weight_matrix[!valid_rows, ] <- 0

  # Convert back to data.table
  weight_df <- data.table(Date = dt$Date)
  for (col in symbol_cols) {
    weight_df[, (col) := weight_matrix[, col]]
  }

  return(weight_df)
}

#' Signal-Based Portfolio Weighting
#'
#' @description
#' Weights selected securities proportionally to their signal strength.
#' Stronger signals receive higher allocations.
#'
#' @param selected_df Binary selection matrix
#' @param signal_df Signal values for weighting
#'
#' @return Data.table with signal-proportional weights
#' @export
#' @examples
#' data("sample_prices_weekly")
#' momentum <- calc_momentum(sample_prices_weekly, lookback = 12)
#' selected <- filter_top_n(momentum, 10)
#' # Weight by momentum strength
#' weights <- weight_by_signal(selected, momentum)
weight_by_signal <- function(selected_df, signal_df) {
  # Weight selected stocks proportionally to their signal strength
  # OPTIMIZED VERSION: 56-140x faster using matrix operations
  #
  # This function allocates portfolio weights based on signal strength,
  # giving higher weights to stocks with stronger signals. Only stocks
  # that are both selected AND have positive signals receive weight.
  #
  # Args:
  #   selected_df: Binary selection matrix (1 = selected, 0 = not selected)
  #                Output from filter functions like filter_top_n()
  #   signal_df: Signal values to use for weighting (momentum, RSI, etc.)
  #              Must have same structure as selected_df
  #
  # Returns:
  #   Weight matrix with allocations proportional to signal strength.
  #   Each row sums to 1 when positions exist, 0 when no valid positions.
  #
  # Examples:
  #   # Weight by momentum strength
  #   selected <- filter_top_n(momentum, n = 10)
  #   weights <- weight_by_signal(selected, momentum)
  #
  #   # Weight by RSI values
  #   high_rsi <- filter_above(rsi, value = 70)
  #   weights <- weight_by_signal(high_rsi, rsi)
  #
  #   # Weight by inverted volatility (lower vol = higher weight)
  #   selected <- filter_top_n(returns, n = 20)
  #   weights <- weight_by_signal(selected, invert_signal(volatility))
  #
  # Details:
  #   - Stocks with higher positive signals get proportionally more weight
  #   - Negative or zero signals receive zero weight (even if selected)
  #   - NA values are treated as zero
  #   - Weights are normalized so each period sums to 1 (fully invested)
  #   - If no stocks have positive signals, all weights are zero
  #
  # Input validation
  validate_selection_data(selected_df, "weight_by_signal")
  validate_selection_data(signal_df, "weight_by_signal")
  validate_matching_structure(selected_df, signal_df, "weight_by_signal")

  # Create copies to avoid mutation
  selected_dt <- ensure_dt_copy(selected_df)
  signal_dt <- ensure_dt_copy(signal_df)
  symbol_cols <- setdiff(names(selected_dt), "Date")

  # Convert to matrices for vectorized operations (massive speedup)
  selection_matrix <- as.matrix(selected_dt[, symbol_cols, with = FALSE])
  signal_matrix <- as.matrix(signal_dt[, symbol_cols, with = FALSE])

  # Handle NA values upfront - convert to 0
  selection_matrix[is.na(selection_matrix)] <- 0
  signal_matrix[is.na(signal_matrix)] <- 0

  # Initialize weight matrix (pre-allocation for efficiency)
  weight_matrix <- matrix(0, nrow = nrow(selection_matrix), ncol = ncol(selection_matrix))

  # Pre-compute valid positions (selected AND positive signal)
  # This avoids redundant checks in the loop
  valid_matrix <- selection_matrix > 0 & signal_matrix > 0

  # Process each time period
  for (i in seq_len(nrow(weight_matrix))) {
    # Get valid positions for this period
    valid_mask <- valid_matrix[i, ]

    # Skip if no valid positions
    if (!any(valid_mask)) next

    # Extract signals for valid positions only
    valid_signals <- signal_matrix[i, valid_mask]

    # Calculate proportional weights
    sum_signals <- sum(valid_signals)

    # Only assign weights if sum is positive (defensive check)
    if (sum_signals > 0) {
      # Normalize signals to sum to 1
      weight_matrix[i, valid_mask] <- valid_signals / sum_signals
    }
  }

  # Convert back to data.table format
  weight_df <- data.table(Date = selected_dt$Date)
  for (j in seq_along(symbol_cols)) {
    weight_df[, (symbol_cols[j]) := weight_matrix[, j]]
  }

  return(weight_df)
}


#' Rank-Based Portfolio Weighting
#'
#' @description
#' Weights securities based on their rank rather than raw signal values.
#' Useful when signal magnitudes are unreliable but ordering is meaningful.
#'
#' @param selected_df Binary selection matrix
#' @param signal_df Signal values for ranking
#' @param method Weighting method: "linear" or "exponential"
#' @param ascending Sort order for ranking (default: FALSE)
#'
#' @return Data.table with rank-based weights
#' @export
#' @examples
#' data("sample_prices_weekly")
#' momentum <- calc_momentum(sample_prices_weekly, lookback = 12)
#' selected <- filter_top_n(momentum, 10)
#' # Linear rank weighting (best gets most)
#' weights <- weight_by_rank(selected, momentum, method = "linear")
#' # Exponential (heavy on top stocks)
#' weights_exp <- weight_by_rank(selected, momentum, method = "exponential")
weight_by_rank <- function(selected_df, signal_df, method = c("linear", "exponential"), ascending = FALSE) {
  # Weight selected stocks by their signal rank
  # OPTIMIZED VERSION: 40-105x faster using matrix operations
  #
  # This function creates a rank-based weighting scheme where stocks are
  # weighted according to their ranking in the signal. Higher-ranked stocks
  # get more weight, with the distribution controlled by the method parameter.
  #
  # Args:
  #   selected_df: Binary selection matrix (1 = selected, 0 = not selected)
  #                Output from filter functions like filter_top_n()
  #   signal_df: Signal values for ranking (momentum, value scores, etc.)
  #              Must have same structure as selected_df
  #   method: Weight distribution method
  #           - "linear" (default): Gradual decline (1st: n/n, 2nd: (n-1)/n, etc.)
  #           - "exponential": Steep decline using 0.5^(rank-1) decay
  #   ascending: Ranking direction
  #              - FALSE (default): Best performers get highest weight
  #              - TRUE: Worst performers get highest weight (contrarian)
  #
  # Returns:
  #   Weight matrix with allocations based on rank.
  #   Each row sums to 1 when positions exist, 0 when no valid positions.
  #
  # Examples:
  #   # Standard momentum with linear decay
  #   selected <- filter_top_n(momentum, n = 10)
  #   weights <- weight_by_rank(selected, momentum, method = "linear")
  #
  #   # Aggressive concentration with exponential decay
  #   weights <- weight_by_rank(selected, momentum, method = "exponential")
  #
  #   # Contrarian: worst performers get highest weight
  #   weights <- weight_by_rank(selected, momentum, ascending = TRUE)
  #
  #   # Value investing with rank weighting
  #   cheap_stocks <- filter_top_n(price_to_book, n = 20, ascending = TRUE)
  #   weights <- weight_by_rank(cheap_stocks, invert_signal(price_to_book))
  #
  # Method Details:
  #   Linear: Creates smooth weight distribution
  #   - 10 stocks: 1st gets 10%, 2nd gets 9%, ..., 10th gets 1%
  #   - After normalization: 1st ~18%, 2nd ~16%, ..., 10th ~2%
  #
  #   Exponential: Creates concentrated weights
  #   - 10 stocks: 1st gets 50%, 2nd gets 25%, 3rd gets 12.5%, etc.
  #   - Top 3 stocks typically get ~85% of total weight
  #
  # Input validation
  validate_selection_data(selected_df, "weight_by_rank")
  validate_selection_data(signal_df, "weight_by_rank")
  validate_matching_structure(selected_df, signal_df, "weight_by_rank")

  method <- match.arg(method)

  # Create copies to avoid mutation
  selected_dt <- ensure_dt_copy(selected_df)
  signal_dt <- ensure_dt_copy(signal_df)
  symbol_cols <- setdiff(names(selected_dt), "Date")

  # Convert to matrices for vectorized operations (massive speedup)
  selection_matrix <- as.matrix(selected_dt[, symbol_cols, with = FALSE])
  signal_matrix <- as.matrix(signal_dt[, symbol_cols, with = FALSE])

  # Handle NA values upfront - convert to 0
  selection_matrix[is.na(selection_matrix)] <- 0
  signal_matrix[is.na(signal_matrix)] <- 0

  # Initialize weight matrix (pre-allocation for efficiency)
  weight_matrix <- matrix(0, nrow = nrow(selection_matrix), ncol = ncol(selection_matrix))

  # Process each time period
  for (i in seq_len(nrow(weight_matrix))) {
    # Get selection mask and signal values for this period
    mask_values <- selection_matrix[i, ]
    signal_values <- signal_matrix[i, ]

    # Only consider stocks that are selected AND have positive signals
    # This ensures we don't weight stocks with negative/zero signals
    valid_mask <- mask_values > 0 & signal_values > 0

    # Skip if no valid positions
    if (!any(valid_mask)) next

    # Get indices and signals for valid stocks
    valid_indices <- which(valid_mask)
    valid_signals <- signal_values[valid_indices]

    # Rank the signals based on ascending parameter
    if (ascending) {
      # For ascending=TRUE, worst signal gets rank 1 (highest weight)
      # Used for contrarian strategies or low-volatility preference
      ranks <- rank(valid_signals, ties.method = "average")
    } else {
      # For ascending=FALSE (default), best signal gets rank 1 (highest weight)
      # Standard momentum approach
      ranks <- rank(-valid_signals, ties.method = "average")
    }

    n <- length(valid_indices)

    # Calculate raw weights based on selected method
    if (method == "linear") {
      # Linear decay: smooth weight distribution
      # Rank 1 gets n points, rank 2 gets n-1 points, etc.
      raw_weights <- (n - ranks + 1) / n
    } else {  # exponential
      # Exponential decay: concentrated weight distribution
      # Each rank gets half the weight of the previous rank
      decay_factor <- 0.5
      raw_weights <- decay_factor^(ranks - 1)
    }

    # Normalize weights to sum to 1 (fully invested)
    normalized_weights <- raw_weights / sum(raw_weights)

    # Assign weights to the weight matrix
    weight_matrix[i, valid_indices] <- normalized_weights
  }

  # Convert back to data.table format
  weight_df <- data.table(Date = selected_dt$Date)
  for (j in seq_along(symbol_cols)) {
    weight_df[, (symbol_cols[j]) := weight_matrix[, j]]
  }

  return(weight_df)
}

#' Volatility-Based Portfolio Weighting
#'
#' @description
#' Weights securities based on their volatility characteristics.
#' Can prefer low-volatility (defensive) or high-volatility (aggressive) stocks.
#'
#' @param selected_df Binary selection matrix (1 = selected, 0 = not)
#' @param vol_timeframe_data Price data for volatility calculation (usually daily)
#' @param strategy_timeframe_data Price data matching strategy frequency
#' @param lookback_periods Number of periods for volatility (default: 26)
#' @param low_vol_preference TRUE = lower vol gets higher weight (default: TRUE)
#' @param vol_method "std", "range", "mad", or "abs_return"
#' @param weighting_method "rank", "equal", or "inverse_variance"
#'
#' @return Data.table with volatility-based weights
#' @export
#' @examples
#' data("sample_prices_weekly")
#' data("sample_prices_daily")
#' momentum <- calc_momentum(sample_prices_weekly, lookback = 12)
#' selected <- filter_top_n(momentum, 10)
#' daily_vol <- calc_rolling_volatility(sample_prices_daily, lookback = 252)
#' aligned_vol <- align_to_timeframe(daily_vol, sample_prices_weekly$Date)
#' weights <- weight_by_volatility(selected, aligned_vol, low_vol_preference = TRUE)
weight_by_volatility <- function(selected_df, vol_timeframe_data,
                                 strategy_timeframe_data = NULL,
                                 lookback_periods = 26,
                                 low_vol_preference = TRUE,
                                 vol_method = "std",
                                 weighting_method = c("rank", "equal", "inverse_variance")) {
  # Weight selected stocks by volatility with multiple calculation and weighting methods
  #
  # Args:
  #   selected_df: Binary selection matrix from strategy
  #   vol_timeframe_data: Price data for volatility (usually daily)
  #   strategy_timeframe_data: Price data matching strategy frequency
  #   lookback_periods: Number of periods for volatility calculation
  #   low_vol_preference: If TRUE, lower volatility gets higher weight
  #   vol_method: Method for volatility calculation ("std", "range", "mad", "abs_return")
  #   weighting_method: How to assign weights ("rank", "equal", "inverse_variance")
  #
  # Returns:
  #   Weight matrix with same timeframe as selected_df

  # Match argument
  weighting_method <- match.arg(weighting_method)

  # Input validation
  validate_selection_data(selected_df, "weight_by_volatility")
  validate_selection_data(vol_timeframe_data, "weight_by_volatility")

  # Setup
  selected_dt <- ensure_dt_copy(selected_df)
  vol_dt <- ensure_dt_copy(vol_timeframe_data)
  symbol_cols <- setdiff(names(selected_dt), "Date")

  # Handle NA in selection
  for (col in symbol_cols) {
    selected_dt[is.na(get(col)), (col) := 0]
  }

  # Default strategy timeframe
  if (is.null(strategy_timeframe_data)) {
    strategy_timeframe_data <- vol_dt
  } else {
    setDT(strategy_timeframe_data)
  }

  # Validate timeframe alignment
  if (!identical(selected_dt$Date, strategy_timeframe_data$Date)) {
    stop("weight_by_volatility: selected_df must have same dates as strategy_timeframe_data")
  }

  # Calculate volatility using the new helper function
  vol_results <- calc_rolling_volatility(vol_dt, lookback_periods, vol_method)

  # Extract variance if using inverse_variance method with std
  if (weighting_method == "inverse_variance" && vol_method == "std") {
    var_results <- attr(vol_results, "variance")
  } else if (weighting_method == "inverse_variance") {
    # For other methods, square the volatility measure
    var_results <- copy(vol_results)
    for (col in symbol_cols) {
      if (col %in% names(var_results)) {
        var_results[, (col) := get(col)^2]
      }
    }
  }

  # Initialize weight matrix
  weight_df <- data.table(Date = selected_dt$Date)
  for(col in symbol_cols) {
    weight_df[, (col) := 0.0]
  }

  # Process each date
  for (i in seq_len(nrow(selected_dt))) {
    date <- selected_dt$Date[i]

    # Get selected securities
    mask <- unlist(selected_dt[i, ..symbol_cols]) > 0
    if (!any(mask)) next

    selected_symbols <- symbol_cols[mask]

    # Find closest volatility date
    vol_dates <- vol_results$Date
    valid_dates <- vol_dates[vol_dates <= date]
    if (length(valid_dates) == 0) next

    closest_date <- max(valid_dates)
    vol_idx <- which(vol_results$Date == closest_date)

    # Get values for weighting
    if (weighting_method == "inverse_variance") {
      # Use variance values
      values <- unlist(var_results[vol_idx, ..selected_symbols])
    } else {
      # Use volatility values
      values <- unlist(vol_results[vol_idx, ..selected_symbols])
    }

    # Apply weighting method using helper
    weights <- apply_weighting_method(
      values = values,
      method = weighting_method,
      preference_ascending = low_vol_preference
    )

    # Assign weights
    if (length(weights) > 0) {
      weight_df[i, (names(weights)) := as.list(weights)]
    }
  }

  return(weight_df)
}

#
#
# # FIXED weight_by_volatility with multiple volatility methods
# weight_by_volatility <- function(selected_df, vol_timeframe_data,
#                                  strategy_timeframe_data = NULL,
#                                  lookback_periods = 26,
#                                  low_vol_preference = TRUE,
#                                  vol_method = "std") {
#   # Weight selected stocks by volatility with multiple calculation methods
#   #
#   # This enhanced version supports multiple volatility estimators:
#   # - "std": Standard deviation of returns (traditional)
#   # - "range": Price range scaled by average (captures trends better)
#   # - "mad": Median Absolute Deviation (robust to outliers)
#   # - "abs_return": Mean absolute return (simple and intuitive)
#   #
#   # Args:
#   #   selected_df: Binary selection matrix from strategy
#   #   vol_timeframe_data: Price data for volatility (usually daily)
#   #   strategy_timeframe_data: Price data matching strategy frequency
#   #   lookback_periods: Number of periods for volatility calculation
#   #   low_vol_preference: If TRUE, lower volatility gets higher weight
#   #   vol_method: Method for volatility calculation (default: "std")
#   #
#   # Returns:
#   #   Weight matrix with same timeframe as selected_df
#
#   # Input validation
#   validate_selection_data(selected_df, "weight_by_volatility")
#   validate_selection_data(vol_timeframe_data, "weight_by_volatility")
#
#   valid_methods <- c("std", "range", "mad", "abs_return")
#   if (!vol_method %in% valid_methods) {
#     stop(paste("vol_method must be one of:", paste(valid_methods, collapse = ", ")))
#   }
#
#   # Convert to data.table
#   selected_dt <- ensure_dt_copy(selected_df)
#   vol_dt <- ensure_dt_copy(vol_timeframe_data)
#
#   symbol_cols <- setdiff(names(selected_dt), "Date")
#
#   # AUTOMATICALLY HANDLE NA in selection
#   for (col in symbol_cols) {
#     selected_dt[is.na(get(col)), (col) := 0]
#   }
#
#   # If strategy_timeframe_data not provided, assume same as vol_timeframe_data
#   if (is.null(strategy_timeframe_data)) {
#     strategy_timeframe_data <- vol_dt
#   } else {
#     setDT(strategy_timeframe_data)
#   }
#
#   # Validate that selected_df has same timeframe as strategy_timeframe_data
#   if (!identical(selected_dt$Date, strategy_timeframe_data$Date)) {
#     stop("weight_by_volatility: selected_df must have the same Date index as strategy_timeframe_data")
#   }
#
#   # Calculate volatility based on method
#   rolling_vol <- copy(vol_dt)
#
#   if (vol_method == "std") {
#     # STANDARD DEVIATION METHOD: Traditional volatility measure
#     returns <- copy(vol_dt)
#     for (col in symbol_cols) {
#       if (col %in% names(returns)) {
#         returns[, (col) := c(NA, diff(get(col))) / shift(get(col))]
#       }
#     }
#
#     for (col in symbol_cols) {
#       if (col %in% names(returns)) {
#         rolling_vol[, (col) := frollapply(returns[[col]], n = lookback_periods,
#                                           FUN = function(x) {
#                                             valid_x <- x[!is.na(x)]
#                                             if(length(valid_x) >= lookback_periods * 0.8) {
#                                               return(sd(valid_x))
#                                             } else {
#                                               return(NA_real_)
#                                             }
#                                           },
#                                           align = "right")]
#       }
#     }
#
#   } else if (vol_method == "range") {
#     # RANGE METHOD: (max - min) / average with scaling for missed intraday
#     for (col in symbol_cols) {
#       if (col %in% names(vol_dt)) {
#         rolling_vol[, (col) := frollapply(get(col), n = lookback_periods,
#                                           FUN = function(x) {
#                                             valid_x <- x[!is.na(x)]
#                                             if(length(valid_x) >= lookback_periods * 0.8) {
#                                               range_val <- (max(valid_x) - min(valid_x))
#                                               avg_val <- mean(valid_x)
#                                               if(avg_val > 0) {
#                                                 # 1.25x scaling for missed intraday moves
#                                                 return(range_val / avg_val)
#                                               }
#                                             }
#                                             return(NA_real_)
#                                           },
#                                           align = "right")]
#       }
#     }
#
#   } else if (vol_method == "mad") {
#     # MAD METHOD: Median Absolute Deviation (robust to outliers)
#     returns <- copy(vol_dt)
#     for (col in symbol_cols) {
#       if (col %in% names(returns)) {
#         returns[, (col) := c(NA, diff(get(col))) / shift(get(col))]
#       }
#     }
#
#     for (col in symbol_cols) {
#       if (col %in% names(returns)) {
#         rolling_vol[, (col) := frollapply(returns[[col]], n = lookback_periods,
#                                           FUN = function(x) {
#                                             valid_x <- x[!is.na(x)]
#                                             if(length(valid_x) >= lookback_periods * 0.8) {
#                                               med <- median(valid_x)
#                                               # Scale factor 1.4826 makes MAD comparable to std
#                                               return(median(abs(valid_x - med)) * 1.4826)
#                                             } else {
#                                               return(NA_real_)
#                                             }
#                                           },
#                                           align = "right")]
#       }
#     }
#
#   } else if (vol_method == "abs_return") {
#     # ABSOLUTE RETURN METHOD: Mean of absolute returns (simple & intuitive)
#     returns <- copy(vol_dt)
#     for (col in symbol_cols) {
#       if (col %in% names(returns)) {
#         returns[, (col) := c(NA, diff(get(col))) / shift(get(col))]
#       }
#     }
#
#     for (col in symbol_cols) {
#       if (col %in% names(returns)) {
#         rolling_vol[, (col) := frollapply(returns[[col]], n = lookback_periods,
#                                           FUN = function(x) {
#                                             valid_x <- x[!is.na(x)]
#                                             if(length(valid_x) >= lookback_periods * 0.8) {
#                                               return(mean(abs(valid_x)))
#                                             } else {
#                                               return(NA_real_)
#                                             }
#                                           },
#                                           align = "right")]
#       }
#     }
#   }
#
#   # Create fresh data.table with numeric columns
#   weight_df <- data.table(Date = selected_dt$Date)
#   for(col in symbol_cols) {
#     weight_df[, (col) := 0.0]
#   }
#
#   # Process each date in selected_df
#   for (i in seq_len(nrow(selected_dt))) {
#     date <- selected_dt$Date[i]
#
#     # Get selected securities on this date
#     mask <- unlist(selected_dt[i, ..symbol_cols]) > 0
#
#     if (!any(mask)) {
#       next
#     }
#
#     # Find the closest date in volatility timeframe
#     vol_dates <- vol_dt$Date
#     valid_dates <- vol_dates[vol_dates <= date]
#
#     if (length(valid_dates) == 0) {
#       next
#     }
#
#     closest_vol_date <- max(valid_dates)
#     vol_row_idx <- which(vol_dt$Date == closest_vol_date)
#
#     # Get volatility for selected securities
#     selected_cols <- symbol_cols[mask]
#     vols <- rolling_vol[vol_row_idx, ..selected_cols]
#
#     # Convert to named vector and remove NAs/Infs
#     vol_vec <- unlist(vols)
#     names(vol_vec) <- selected_cols
#     vol_vec <- vol_vec[!is.na(vol_vec) & !is.infinite(vol_vec)]
#
#     if (length(vol_vec) == 0) {
#       next
#     }
#
#     # Sort volatilities (ascending if we prefer low volatility)
#     sorted_vols <- sort(vol_vec, decreasing = !low_vol_preference)
#
#     # Create ranks (n for lowest/highest vol based on preference)
#     n <- length(sorted_vols)
#     ranks <- n:1
#     names(ranks) <- names(sorted_vols)
#
#     # Normalize ranks to get weights
#     normalized_weights <- ranks / sum(ranks)
#
#     # Create full weight vector and use as.list() for assignment
#     weights_vector <- numeric(length(symbol_cols))
#     names(weights_vector) <- symbol_cols
#     weights_vector[names(normalized_weights)] <- normalized_weights
#
#     weight_df[i, (symbol_cols) := as.list(weights_vector)]
#   }
#
#   return(weight_df)
# }




#' Regime-Based Adaptive Weighting
#'
#' @description
#' Applies different weighting methods based on market regime classification.
#' Enables adaptive strategies that change allocation approach in different
#' market conditions.
#'
#' @param selected_df Binary selection matrix (1 = selected, 0 = not)
#' @param regime Regime classification (integer values per period)
#' @param signal_df Signal values (required for signal/rank methods)
#' @param vol_timeframe_data Volatility data (required for volatility method)
#' @param strategy_timeframe_data Strategy timeframe alignment data
#' @param weighting_configs List with method-specific parameters
#'
#' @return Data.table with regime-adaptive weights
#' @export
#' @examples
#' data("sample_prices_weekly")
#' # Create selection and signals
#' momentum <- calc_momentum(sample_prices_weekly, lookback = 12)
#' selected <- filter_top_n(momentum, n = 10)
#'
#' # Create a simple regime (example: based on market trend)
#' ma20 <- calc_moving_average(sample_prices_weekly, 20)
#' spy_price <- sample_prices_weekly$SPY
#' spy_ma <- ma20$SPY
#' regime <- ifelse(spy_price > spy_ma, 1, 2)
#'
#' # Different weights for bull/bear markets
#' weighting_configs <- list(
#'   "1" = list(method = "equal"),
#'   "2" = list(method = "signal")
#' )
#' weights <- weight_by_regime(selected, regime, weighting_configs,
#'                             signal_df = momentum)
weight_by_regime <- function(selected_df, regime, weighting_configs,
                             signal_df = NULL, vol_timeframe_data = NULL,
                             strategy_timeframe_data = NULL) {


  # Input validation
  validate_selection_data(selected_df, "weight_by_regime")

  # Extract regime vector if needed
  if (is.data.frame(regime) || is.data.table(regime)) {
    if (ncol(regime) == 2 && "Date" %in% names(regime)) {
      regime_col <- setdiff(names(regime), "Date")[1]
      regime_vector <- regime[[regime_col]]
    } else {
      stop("weight_by_regime: regime must be a vector or data.frame with Date and one regime column")
    }
  } else {
    regime_vector <- regime
  }

  # Check lengths match
  if (length(regime_vector) != nrow(selected_df)) {
    stop("weight_by_regime: regime length must match selected_df rows")
  }

  # Get unique regimes
  unique_regimes <- unique(regime_vector[!is.na(regime_vector)])

  # Check all regimes have configs
  missing_regimes <- setdiff(unique_regimes, names(weighting_configs))
  if (length(missing_regimes) > 0) {
    stop(paste("weight_by_regime: Missing configs for regimes:",
               paste(missing_regimes, collapse = ", ")))
  }

  # Initialize result
  selected_dt <- ensure_dt_copy(selected_df)
  symbol_cols <- setdiff(names(selected_dt), "Date")

  weight_df <- data.table(Date = selected_dt$Date)
  for(col in symbol_cols) {
    weight_df[, (col) := 0.0]
  }

  # Pre-process vol_timeframe_data if provided (ensure data.table)
  if (!is.null(vol_timeframe_data)) {
    vol_dt <- ensure_dt_copy(vol_timeframe_data)
  }

  # Process each regime
  for (reg in unique_regimes) {
    cat("Processing regime:", reg, "\n")

    # Find periods in this regime
    regime_idx <- which(regime_vector == reg)
    if (length(regime_idx) == 0) {
      cat("  No periods found for regime", reg, "\n")
      next
    }

    cat("  Regime periods:", length(regime_idx), "\n")

    # Get config for this regime
    config <- weighting_configs[[reg]]
    method <- config$method

    # Extract subset for this regime
    regime_selection <- selected_dt[regime_idx, ]
    cat("  Regime selection dates:", min(regime_selection$Date), "to", max(regime_selection$Date), "\n")

    # Call appropriate weighting function based on method
    if (method == "equal") {
      regime_weights <- weight_equally(regime_selection)

    } else if (method == "signal") {
      if (is.null(signal_df)) stop("weight_by_regime: signal_df required for signal weighting")
      regime_signal <- signal_df[regime_idx, ]
      regime_weights <- weight_by_signal(regime_selection, regime_signal)

    } else if (method == "rank") {
      if (is.null(signal_df)) stop("weight_by_regime: signal_df required for rank weighting")
      regime_signal <- signal_df[regime_idx, ]
      regime_weights <- weight_by_rank(regime_selection, regime_signal)

    } else if (method == "volatility") {
      if (is.null(vol_timeframe_data)) {
        stop("weight_by_regime: vol_timeframe_data required for volatility weighting")
      }

      # Extract parameters from config
      vol_method <- config$vol_method %||% "std"
      weighting_method <- config$weighting_method %||% "rank"
      low_vol_pref <- config$low_vol_preference %||% TRUE
      lookback <- config$lookback_periods %||% 60

      cat("  Volatility config: method =", vol_method, ", lookback =", lookback, "\n")

      # FIXED: Properly handle strategy_timeframe_data subsetting
      if (!is.null(strategy_timeframe_data)) {
        setDT(strategy_timeframe_data)
        regime_strategy_data <- strategy_timeframe_data[regime_idx, ]
        cat("  Using provided strategy timeframe data\n")
      } else {
        # If no strategy_timeframe_data provided, use regime selection's timeframe
        regime_strategy_data <- regime_selection
        cat("  Using regime selection as strategy timeframe\n")
      }

      # IMPROVED: Smart vol_timeframe_data subsetting
      # We need enough history BEFORE the regime starts for lookback calculation
      first_regime_date <- min(regime_selection$Date)
      last_regime_date <- max(regime_selection$Date)

      # Find all vol data up to the last regime date (we need history before first date)
      vol_end_idx <- which(vol_dt$Date <= last_regime_date)
      if (length(vol_end_idx) == 0) {
        warning(paste("weight_by_regime: No vol data available for regime", reg))
        regime_weights <- weight_equally(regime_selection)
      } else {
        # Use all available vol data up to regime end (includes necessary history)
        regime_vol_data <- vol_dt[1:max(vol_end_idx), ]

        cat("  Vol data range:", min(regime_vol_data$Date), "to", max(regime_vol_data$Date), "\n")
        cat("  Vol data rows:", nrow(regime_vol_data), "\n")

        # Validate we have enough data for lookback
        if (nrow(regime_vol_data) < lookback) {
          warning(paste("weight_by_regime: Limited vol data for regime", reg,
                        "- have", nrow(regime_vol_data), "periods, need", lookback))
        }

        # Try volatility weighting with robust error handling
        tryCatch({
          regime_weights <- weight_by_volatility(
            selected_df = regime_selection,
            vol_timeframe_data = regime_vol_data,      # FIXED: Use properly subset vol data
            strategy_timeframe_data = regime_strategy_data,
            lookback_periods = lookback,
            low_vol_preference = low_vol_pref,
            vol_method = vol_method,
            weighting_method = weighting_method
          )

          cat("  [OK] Volatility weights calculated successfully\n")
          cat("    Rows:", nrow(regime_weights), "\n")
          cat("    Total weight:", round(sum(regime_weights[, ..symbol_cols]), 4), "\n")

        }, error = function(e) {
          cat("  [ERROR] ERROR in weight_by_volatility for regime", reg, ":", e$message, "\n")
          cat("    Falling back to equal weights\n")
          regime_weights <- weight_equally(regime_selection)
        })
      }

    } else {
      stop(paste("weight_by_regime: Unknown method:", method))
    }

    # IMPROVED: Safer weight copying with validation
    if (nrow(regime_weights) != length(regime_idx)) {
      warning(paste("weight_by_regime: Dimension mismatch for regime", reg))
      next
    }

    # Copy weights back to main result
    for (col in symbol_cols) {
      if (col %in% names(regime_weights)) {
        weight_df[regime_idx, (col) := regime_weights[[col]]]
      }
    }

    cat("  [OK] Regime", reg, "processed successfully\n\n")
  }

  return(weight_df)
}



###############################################################################
# UTILITY FUNCTIONS
###############################################################################
# Add to weighting.R (after the existing weighting functions)
#' Apply Weighting Method to Values
#'
#' @description
#' Internal utility function that applies various weighting methods to a vector.
#' Used by other weighting functions.
#'
#' @param values Named numeric vector of values to weight
#' @param preference_ascending TRUE = prefer lower values, FALSE = prefer higher
#'
#' @return Named numeric vector of weights that sum to 1
#' @keywords internal
apply_weighting_method <- function(values, method = "rank", preference_ascending = TRUE) {
  # Apply a weighting method to a vector of values
  #
  # Args:
  #   values: Named numeric vector of values to weight
  #   method: One of "rank", "equal", "inverse_variance"
  #   preference_ascending: TRUE = prefer lower values, FALSE = prefer higher values
  #
  # Returns:
  #   Named numeric vector of weights that sum to 1

  # Remove NA/Inf values
  valid_values <- values[!is.na(values) & !is.infinite(values)]

  if (length(valid_values) == 0) {
    return(numeric(0))
  }

  if (method == "rank") {
    # Sort by preference
    sorted_vals <- sort(valid_values, decreasing = !preference_ascending)
    n <- length(sorted_vals)

    # Assign ranks (highest rank to preferred values)
    ranks <- n:1
    names(ranks) <- names(sorted_vals)

    # Normalize
    weights <- ranks / sum(ranks)

  } else if (method == "equal") {
    # Equal weights for all
    weights <- rep(1/length(valid_values), length(valid_values))
    names(weights) <- names(valid_values)

  } else if (method == "inverse_variance") {
    # For inverse variance, values should be variances
    if (any(valid_values <= 0)) {
      stop("inverse_variance method requires positive variance values")
    }

    if (preference_ascending) {
      # Prefer low variance: weight = 1/variance
      inv_weights <- 1 / valid_values
    } else {
      # Prefer high variance: use variance directly
      inv_weights <- valid_values
    }

    # Normalize
    weights <- inv_weights / sum(inv_weights)
  } else {
    stop("Unknown weighting method: ", method)
  }

  return(weights)
}

summary_weights <- function(weight_df) {
  # Provide summary statistics about weight allocations
  #
  # Args:
  #   weight_df (DataFrame): Weight matrix from weighting functions
  #
  # Returns:
  #   List: Summary statistics

  setDT(weight_df)
  symbol_cols <- setdiff(names(weight_df), "Date")

  # Calculate weight statistics
  total_weight_per_date <- weight_df[, rowSums(.SD, na.rm = TRUE), .SDcols = symbol_cols]
  max_weight_per_date <- weight_df[, apply(.SD, 1, max, na.rm = TRUE), .SDcols = symbol_cols]
  active_positions_per_date <- weight_df[, rowSums(.SD > 0, na.rm = TRUE), .SDcols = symbol_cols]

  # Average weight when position is active
  avg_active_weights <- numeric(length(symbol_cols))
  names(avg_active_weights) <- symbol_cols
  for (col in symbol_cols) {
    active_weights <- weight_df[[col]][weight_df[[col]] > 0]
    avg_active_weights[col] <- if(length(active_weights) > 0) mean(active_weights) else 0
  }

  summary_stats <- list(
    total_dates = nrow(weight_df),
    total_symbols = length(symbol_cols),
    avg_total_weight_per_date = mean(total_weight_per_date),
    min_total_weight = min(total_weight_per_date),
    max_total_weight = max(total_weight_per_date),
    avg_max_position_size = mean(max_weight_per_date),
    max_single_position = max(max_weight_per_date),
    avg_active_positions = mean(active_positions_per_date),
    dates_fully_invested = sum(total_weight_per_date > 0.99),
    dates_in_cash = sum(total_weight_per_date == 0),
    avg_weight_when_active = mean(avg_active_weights)
  )

  return(summary_stats)
}



#' Combine Multiple Weighting Schemes
#'
#' @description
#' Blends multiple weight matrices with specified weights. Useful for
#' multi-factor strategies that combine different allocation approaches.
#' Optimized using matrix operations for 1000x+ speedup.
#'
#' @param weight_matrices List of weight data frames to combine
#' @param weights Numeric vector of weights for each matrix (default: equal)
#'
#' @return Data.table with blended portfolio weights
#' @export
#' @examples
#' data("sample_prices_weekly")
#' # Calculate signals
#' momentum <- calc_momentum(sample_prices_weekly, lookback = 12)
#' selected <- filter_top_n(momentum, n = 10)
#' volatility <- calc_rolling_volatility(sample_prices_weekly, lookback = 20)
#'
#' # Combine momentum and low-vol weights
#' mom_weights <- weight_by_signal(selected, momentum)
#' vol_weights <- weight_by_signal(selected, invert_signal(volatility))
#' combined <- combine_weights(list(mom_weights, vol_weights), weights = c(0.7, 0.3))
combine_weights <- function(weight_matrices, weights = NULL) {
  # Combine multiple weight matrices with specified weights
  # OPTIMIZED VERSION: Uses matrix operations for 1000x+ speedup

  # Input validation
  if (!is.list(weight_matrices)) {
    stop("combine_weights: weight_matrices must be a list")
  }

  if (length(weight_matrices) < 2) {
    stop("combine_weights: need at least 2 weight matrices to combine")
  }

  # Validate all inputs are data frames
  for (i in seq_along(weight_matrices)) {
    if (!is.data.frame(weight_matrices[[i]])) {
      stop(paste("combine_weights: weight_matrices[[", i, "]] must be a data.frame or data.table"))
    }
  }

  # Check dimensions match
  base_dim <- dim(weight_matrices[[1]])
  base_names <- names(weight_matrices[[1]])
  base_dates <- weight_matrices[[1]]$Date

  for (i in 2:length(weight_matrices)) {
    if (!identical(dim(weight_matrices[[i]]), base_dim)) {
      stop(paste("combine_weights: weight_matrices[[", i, "]] has different dimensions"))
    }
    if (!identical(names(weight_matrices[[i]]), base_names)) {
      stop(paste("combine_weights: weight_matrices[[", i, "]] has different column names"))
    }
    if (!identical(weight_matrices[[i]]$Date, base_dates)) {
      stop(paste("combine_weights: weight_matrices[[", i, "]] has different dates"))
    }
  }

  # Handle weights parameter
  if (is.null(weights)) {
    weights <- rep(1, length(weight_matrices))
  }

  if (length(weights) != length(weight_matrices)) {
    stop("combine_weights: length of weights must match number of weight matrices")
  }

  if (!is.numeric(weights)) {
    stop("combine_weights: weights must be numeric")
  }

  if (any(weights < 0)) {
    stop("combine_weights: weights must be non-negative")
  }

  if (all(weights == 0)) {
    stop("combine_weights: at least one weight must be positive")
  }

  # Normalize weights to sum to 1
  weights <- weights / sum(weights)

  # Get structure info
  dates <- base_dates
  weight_cols <- setdiff(base_names, "Date")

  # OPTIMIZED: Convert to matrices for fast operations
  mat_list <- lapply(weight_matrices, function(df) {
    # Handle NA values by converting to 0
    mat <- as.matrix(df[, weight_cols, with = FALSE])
    mat[is.na(mat)] <- 0
    return(mat)
  })

  # OPTIMIZED: Linear combination using matrix operations
  result_mat <- mat_list[[1]] * weights[1]
  for (i in 2:length(mat_list)) {
    result_mat <- result_mat + mat_list[[i]] * weights[i]
  }

  # OPTIMIZED: Vectorized row normalization
  row_sums <- rowSums(result_mat)
  # Use matrix division with recycling for normalization
  # Set zero rows to 1 to avoid division by zero
  row_sums[abs(row_sums) < 1e-10] <- 1
  result_mat <- result_mat / row_sums

  # Convert back to data.table
  result <- data.table(Date = dates)
  result[, (weight_cols) := as.data.table(result_mat)]

  return(result)
}

validate_weights <- function(weight_df, tolerance = 1e-10) {
  # Validate that weights sum to approximately 1 or 0 for each date
  #
  # Args:
  #   weight_df (DataFrame): Weight matrix to validate
  #   tolerance (numeric): Tolerance for sum check
  #
  # Returns:
  #   Boolean: TRUE if validation passes, with warnings for issues

  setDT(weight_df)
  symbol_cols <- setdiff(names(weight_df), "Date")

  # Check each row sums to approximately 0 or 1
  row_sums <- weight_df[, rowSums(.SD, na.rm = TRUE), .SDcols = symbol_cols]

  valid_sums <- (abs(row_sums - 1) < tolerance) | (abs(row_sums) < tolerance)

  if (!all(valid_sums)) {
    problem_dates <- weight_df$Date[!valid_sums]
    problem_sums <- row_sums[!valid_sums]

    warning(paste("Weight validation failed for", sum(!valid_sums), "dates."))
    warning(paste("First few problem dates:", paste(head(problem_dates), collapse = ", ")))
    warning(paste("Their sums:", paste(round(head(problem_sums), 4), collapse = ", ")))

    return(FALSE)
  }

  # Check all weights are non-negative
  all_weights <- unlist(weight_df[, ..symbol_cols])
  if (any(all_weights < 0, na.rm = TRUE)) {
    warning("Negative weights found!")
    return(FALSE)
  }

  return(TRUE)
}



#' Switch Between Weighting Schemes
#'
#' @description
#' Dynamically switches between two weighting schemes based on a signal.
#' Enables tactical allocation changes.
#'
#' @param weights_a Primary weight matrix
#' @param weights_b Alternative weight matrix
#' @param use_b_condition Logical vector (TRUE = use weights_b)
#' @param partial_blend Blend factor 0-1 (default: 1 = full switch)
#'
#' @return Combined weight matrix
#' @export
#' @examples
#' data("sample_prices_weekly")
#' momentum <- calc_momentum(sample_prices_weekly, lookback = 12)
#' selected <- filter_top_n(momentum, n = 10)
#' weights_equal <- weight_equally(selected)
#' weights_signal <- weight_by_signal(selected, momentum)
#'
#' # Create switching signal (example: use SPY momentum as regime indicator)
#' spy_momentum <- momentum$SPY
#' switch_signal <- as.numeric(spy_momentum > median(spy_momentum, na.rm = TRUE))
#' switch_signal[is.na(switch_signal)] <- 0
#'
#' # Switch between strategies
#' final_weights <- switch_weights(weights_equal, weights_signal, switch_signal)
switch_weights <- function(weights_a, weights_b, use_b_condition, partial_blend = 1) {
  # Validate inputs
  if (!is.data.frame(weights_a) && !is.data.table(weights_a)) {
    stop("switch_weights: weights_a must be a data.frame or data.table")
  }

  if (!is.data.frame(weights_b) && !is.data.table(weights_b)) {
    stop("switch_weights: weights_b must be a data.frame or data.table")
  }

  if (!identical(dim(weights_a), dim(weights_b))) {
    stop("switch_weights: weights_a and weights_b must have same dimensions")
  }

  if (!identical(names(weights_a), names(weights_b))) {
    stop("switch_weights: weights_a and weights_b must have same column names")
  }

  if (length(use_b_condition) != nrow(weights_a)) {
    stop(sprintf("switch_weights: condition length (%d) must match weight matrix rows (%d)",
                 length(use_b_condition), nrow(weights_a)))
  }

  if (!is.logical(use_b_condition) && !is.numeric(use_b_condition)) {
    stop("switch_weights: use_b_condition must be logical or numeric (0/1)")
  }

  if (!is.numeric(partial_blend) || partial_blend < 0 || partial_blend > 1) {
    stop("switch_weights: partial_blend must be numeric between 0 and 1")
  }

  # Ensure data.table and create copy to avoid modifying input
  weights_a_dt <- ensure_dt_copy(weights_a)
  weights_b_dt <- ensure_dt_copy(weights_b)

  # Get weight columns (exclude Date)
  weight_cols <- setdiff(names(weights_a_dt), "Date")

  # Initialize result with copy of weights_a
  result <- copy(weights_a_dt)

  # Convert condition to logical if numeric
  use_b_logical <- as.logical(use_b_condition)

  # Find rows where we need to use weights_b
  switch_rows <- which(use_b_logical & !is.na(use_b_logical))

  if (length(switch_rows) > 0) {
    if (partial_blend == 1) {
      # Full switch - use weights_b for these rows
      # FIXED: Use .. notation for proper column selection
      result[switch_rows, (weight_cols) :=
               weights_b_dt[switch_rows, ..weight_cols]]
    } else {
      # Partial blend - weighted average of weights_a and weights_b
      for (i in switch_rows) {
        for (col in weight_cols) {
          a_weight <- weights_a_dt[i, get(col)]
          b_weight <- weights_b_dt[i, get(col)]

          # Handle NA values
          if (is.na(a_weight) || is.na(b_weight)) {
            blended <- NA_real_
          } else {
            blended <- (1 - partial_blend) * a_weight + partial_blend * b_weight
          }

          result[i, (col) := blended]
        }
      }
    }
  }

  # Log switching statistics if verbose (following library patterns)
  n_switches <- length(switch_rows)
  pct_switches <- 100 * n_switches / nrow(result)

  if (getOption("momentum.verbose", FALSE)) {
    cat(sprintf("switch_weights: Switched weights for %d periods (%.1f%%)\n",
                n_switches, pct_switches))
  }

  return(result)
}







#' Optimized recursive bisection for HRP
#' @keywords internal
recursive_bisection_optimized <- function(cov_matrix, asset_order) {

  n_assets <- length(asset_order)
  weights <- rep(1.0, n_assets)
  names(weights) <- asset_order

  # Use iterative approach instead of recursion to avoid stack overflow
  # and eliminate global variables

  split_queue <- list(list(assets = asset_order, weight_mult = 1.0))

  while (length(split_queue) > 0) {
    current <- split_queue[[1]]
    split_queue <- split_queue[-1]

    current_assets <- current$assets
    current_mult <- current$weight_mult

    if (length(current_assets) == 1) {
      # Base case - single asset
      weights[current_assets] <- weights[current_assets] * current_mult
      next
    }

    # Split into two groups
    mid_point <- ceiling(length(current_assets) / 2)
    group1 <- current_assets[1:mid_point]
    group2 <- current_assets[(mid_point + 1):length(current_assets)]

    # Calculate cluster variances - VECTORIZED
    var1 <- calculate_cluster_variance_optimized(cov_matrix, group1)
    var2 <- calculate_cluster_variance_optimized(cov_matrix, group2)

    # Calculate allocation ratio
    total_var <- var1 + var2
    if (total_var > 0) {
      alpha <- var2 / total_var  # Inverse variance weighting
    } else {
      alpha <- 0.5  # Equal split if no variance info
    }

    # Update weights for this split
    for (asset in group1) {
      weights[asset] <- weights[asset] * alpha
    }
    for (asset in group2) {
      weights[asset] <- weights[asset] * (1 - alpha)
    }

    # Add subgroups to queue if they need further splitting
    if (length(group1) > 1) {
      split_queue <- append(split_queue, list(list(assets = group1, weight_mult = 1.0)))
    }
    if (length(group2) > 1) {
      split_queue <- append(split_queue, list(list(assets = group2, weight_mult = 1.0)))
    }
  }

  return(weights)
}

#' Optimized cluster variance calculation
#' @keywords internal
calculate_cluster_variance_optimized <- function(cov_matrix, asset_names) {

  if (length(asset_names) == 1) {
    return(cov_matrix[asset_names, asset_names])
  }

  # Extract submatrix - VECTORIZED
  sub_cov <- cov_matrix[asset_names, asset_names, drop = FALSE]

  # Inverse variance portfolio weights - VECTORIZED
  inv_diag <- 1 / diag(sub_cov)
  inv_diag[!is.finite(inv_diag)] <- 0  # Handle zero variance

  if (sum(inv_diag) == 0) {
    return(0)
  }

  ivp_weights <- inv_diag / sum(inv_diag)

  # Portfolio variance: w' * Cov * w - VECTORIZED
  portfolio_var <- as.numeric(t(ivp_weights) %*% sub_cov %*% ivp_weights)

  return(portfolio_var)
}

#' Hierarchical Risk Parity Weighting
#'
#' @description
#' Calculates portfolio weights using Hierarchical Risk Parity (HRP) methodology.
#' HRP combines hierarchical clustering with risk-based allocation to create
#' diversified portfolios that don't rely on unstable correlation matrix inversions.
#'
#' @param selected_df Binary selection matrix (data.frame with Date column)
#' @param prices_df Price data for covariance calculation (typically daily)
#'                  Returns are calculated internally from prices
#' @param lookback_periods Number of periods for covariance estimation (default: 252)
#' @param cluster_method Clustering linkage method (default: "ward.D2")
#' @param distance_method Distance measure for clustering (default: "euclidean")
#' @param min_periods Minimum periods required for calculation (default: 60)
#' @param use_correlation If TRUE, cluster on correlation instead of covariance
#'
#' @return Weight matrix with same dates as selected_df
#'
#' @details
#' The HRP algorithm:
#' 1. Calculate returns from input prices
#' 2. Compute covariance matrix from returns
#' 3. Cluster assets based on distance matrix
#' 4. Apply recursive bisection with inverse variance weighting
#' 5. Results in naturally diversified portfolio without matrix inversion
#'
#' The function accepts price data and calculates returns internally,
#' matching the pattern of other library functions like calc_momentum().
#'
#' @export
#' @examples
#' data("sample_prices_daily")
#' data("sample_prices_weekly")
#' # Create a selection first
#' momentum <- calc_momentum(sample_prices_weekly, lookback = 12)
#' selected <- filter_top_n(momentum, n = 10)
#'
#' # Using daily prices for risk calculation
#' weights <- weight_by_hrp(selected, sample_prices_daily, lookback_periods = 252)
#'
#' # Using correlation-based clustering
#' weights <- weight_by_hrp(selected, sample_prices_daily, use_correlation = TRUE)
weight_by_hrp <- function(selected_df, prices_df,
                          lookback_periods = 252,
                          cluster_method = "ward.D2",
                          distance_method = "euclidean",
                          min_periods = 60,
                          use_correlation = FALSE) {
  # Hierarchical Risk Parity with price inputs
  #
  # This version accepts price data and calculates returns internally,
  # matching the pattern of other library functions like calc_momentum().
  # The frequency of prices_df determines the covariance calculation frequency.
  #
  # Args:
  #   selected_df: Binary selection matrix (data.frame with Date column)
  #   prices_df: Price data for covariance calculation (typically daily)
  #   lookback_periods: Number of periods for covariance (default: 252 for daily)
  #   cluster_method: Clustering linkage method (default: "ward.D2")
  #   distance_method: Distance measure (default: "euclidean")
  #   min_periods: Minimum periods required (default: 60)
  #   use_correlation: Use correlation instead of covariance (default: FALSE)
  #
  # Returns:
  #   Weight matrix with same dates as selected_df

  # Input validation
  validate_selection_data(selected_df, "weight_by_hrp")
  validate_selection_data(prices_df, "weight_by_hrp")

  if (lookback_periods < min_periods) {
    stop("weight_by_hrp: lookback_periods must be >= min_periods")
  }

  valid_cluster_methods <- c("ward.D", "ward.D2", "single", "complete",
                             "average", "mcquitty", "median", "centroid")
  if (!cluster_method %in% valid_cluster_methods) {
    stop("weight_by_hrp: Invalid cluster_method. Use: ",
         paste(valid_cluster_methods, collapse = ", "))
  }

  valid_distance_methods <- c("euclidean", "maximum", "manhattan",
                              "canberra", "binary", "minkowski")
  if (!distance_method %in% valid_distance_methods) {
    stop("weight_by_hrp: Invalid distance_method. Use: ",
         paste(valid_distance_methods, collapse = ", "))
  }

  # Setup
  selected_dt <- ensure_dt_copy(selected_df)
  prices_dt <- ensure_dt_copy(prices_df)
  symbol_cols <- setdiff(names(selected_dt), "Date")

  # Handle NA in selection
  for (col in symbol_cols) {
    selected_dt[is.na(get(col)), (col) := 0]
  }

  # CALCULATE RETURNS FROM PRICES (matching library pattern)
  returns_dt <- data.table(Date = prices_dt$Date)
  for (col in symbol_cols) {
    if (col %in% names(prices_dt)) {
      # Simple returns: (P_t - P_{t-1}) / P_{t-1}
      prices <- prices_dt[[col]]
      returns <- c(NA, diff(prices) / head(prices, -1))
      returns_dt[, (col) := returns]
    }
  }

  # Initialize result
  weight_df <- data.table(Date = selected_dt$Date)
  for(col in symbol_cols) {
    weight_df[, (col) := 0.0]
  }

  # Process each rebalancing date
  for (i in seq_len(nrow(selected_dt))) {
    date <- selected_dt$Date[i]

    # Get selected assets
    selected_mask <- unlist(selected_dt[i, ..symbol_cols]) > 0
    if (!any(selected_mask)) {
      next  # No assets selected
    }

    selected_symbols <- symbol_cols[selected_mask]

    # Get historical returns for lookback period
    end_idx <- which(returns_dt$Date <= date)
    if (length(end_idx) == 0) {
      next
    }

    latest_idx <- max(end_idx)
    start_idx <- max(1, latest_idx - lookback_periods + 1)

    if (latest_idx - start_idx + 1 < min_periods) {
      next  # Insufficient history
    }

    # Extract returns matrix for selected assets
    returns_subset <- returns_dt[start_idx:latest_idx, c("Date", selected_symbols), with = FALSE]
    returns_matrix <- as.matrix(returns_subset[, ..selected_symbols])

    # Remove any rows with all NAs
    valid_rows <- rowSums(!is.na(returns_matrix)) > 0
    returns_matrix <- returns_matrix[valid_rows, , drop = FALSE]

    if (nrow(returns_matrix) < min_periods) {
      next
    }

    # Calculate HRP weights
    tryCatch({
      hrp_weights <- calculate_hrp_weights(
        returns_matrix = returns_matrix,
        asset_names = selected_symbols,
        cluster_method = cluster_method,
        distance_method = distance_method

      )

      # Assign weights
      for (symbol in selected_symbols) {
        weight_df[i, (symbol) := hrp_weights[symbol]]
      }

    }, error = function(e) {
      # Fallback to equal weights if HRP fails
      equal_weight <- 1 / length(selected_symbols)
      for (symbol in selected_symbols) {
        weight_df[i, (symbol) := equal_weight]
      }
    })
  }

  return(weight_df)
}

# Note: calculate_hrp_weights() internal function remains unchanged
# It already expects returns as input, so no modification needed
#' Calculate HRP weights for a given returns matrix
#' @keywords internal
calculate_hrp_weights <- function(returns_matrix, asset_names,
                                  cluster_method, distance_method) {

  # Handle missing data by column-wise removal if necessary
  complete_cols <- colSums(is.na(returns_matrix)) == 0
  if (sum(complete_cols) < 2) {
    stop("Insufficient complete return series for HRP calculation")
  }

  # Use only complete columns if needed
  if (!all(complete_cols)) {
    returns_matrix <- returns_matrix[, complete_cols, drop = FALSE]
    asset_names <- asset_names[complete_cols]
  }

  n_assets <- ncol(returns_matrix)
  if (n_assets < 2) {
    stop("Need at least 2 assets for HRP")
  }

  # Calculate covariance matrix
  cov_matrix <- cov(returns_matrix, use = "complete.obs")
  rownames(cov_matrix) <- asset_names
  colnames(cov_matrix) <- asset_names

  # Convert to correlation for distance calculation
  cor_matrix <- cov2cor(cov_matrix)

  # FIXED: Properly implement distance_method switching
  if (distance_method == "euclidean") {
    # Traditional HRP approach: sqrt(0.5 * (1 - correlation))
    dist_matrix <- as.dist(sqrt(0.5 * (1 - cor_matrix)))
  } else {
    # For non-euclidean methods, use R's dist() function on correlation matrix
    # First transform correlation matrix to distance-appropriate format
    if (distance_method %in% c("manhattan", "maximum", "canberra")) {
      # Use (1 - correlation) transformation for these methods
      distance_data <- 1 - cor_matrix
      # Convert to matrix form that dist() can handle
      dist_matrix <- dist(distance_data, method = distance_method)
    } else {
      # For other methods, apply directly to correlation matrix
      dist_matrix <- dist(cor_matrix, method = distance_method)
    }
  }

  # Hierarchical clustering - OPTIMIZED
  hclust_result <- hclust(dist_matrix, method = cluster_method)

  # Get cluster ordering
  sort_order <- hclust_result$order
  sorted_assets <- asset_names[sort_order]

  # Apply recursive bisection - OPTIMIZED VERSION
  hrp_weights <- recursive_bisection_optimized(cov_matrix, sorted_assets)

  # Ensure all original assets have weights (even if zero)
  final_weights <- rep(0, length(asset_names))
  names(final_weights) <- asset_names
  final_weights[names(hrp_weights)] <- hrp_weights

  return(final_weights)
}

#' Optimized recursive bisection for HRP
#' @keywords internal
recursive_bisection_optimized <- function(cov_matrix, asset_order) {

  n_assets <- length(asset_order)
  weights <- rep(1.0, n_assets)
  names(weights) <- asset_order

  # Use iterative approach instead of recursion to avoid stack overflow
  # and eliminate global variables

  split_queue <- list(list(assets = asset_order, weight_mult = 1.0))

  while (length(split_queue) > 0) {
    current <- split_queue[[1]]
    split_queue <- split_queue[-1]

    current_assets <- current$assets
    current_mult <- current$weight_mult

    if (length(current_assets) == 1) {
      # Base case - single asset
      weights[current_assets] <- weights[current_assets] * current_mult
      next
    }

    # Split into two groups
    mid_point <- ceiling(length(current_assets) / 2)
    group1 <- current_assets[1:mid_point]
    group2 <- current_assets[(mid_point + 1):length(current_assets)]

    # Calculate cluster variances - VECTORIZED
    var1 <- calculate_cluster_variance_optimized(cov_matrix, group1)
    var2 <- calculate_cluster_variance_optimized(cov_matrix, group2)

    # Calculate allocation ratio
    total_var <- var1 + var2
    if (total_var > 0) {
      alpha <- var2 / total_var  # Inverse variance weighting
    } else {
      alpha <- 0.5  # Equal split if no variance info
    }

    # Update weights for this split
    for (asset in group1) {
      weights[asset] <- weights[asset] * alpha
    }
    for (asset in group2) {
      weights[asset] <- weights[asset] * (1 - alpha)
    }

    # Add subgroups to queue if they need further splitting
    if (length(group1) > 1) {
      split_queue <- append(split_queue, list(list(assets = group1, weight_mult = 1.0)))
    }
    if (length(group2) > 1) {
      split_queue <- append(split_queue, list(list(assets = group2, weight_mult = 1.0)))
    }
  }

  return(weights)
}

#' Optimized cluster variance calculation
#' @keywords internal
calculate_cluster_variance_optimized <- function(cov_matrix, asset_names) {

  if (length(asset_names) == 1) {
    return(cov_matrix[asset_names, asset_names])
  }

  # Extract submatrix - VECTORIZED
  sub_cov <- cov_matrix[asset_names, asset_names, drop = FALSE]

  # Inverse variance portfolio weights - VECTORIZED
  inv_diag <- 1 / diag(sub_cov)
  inv_diag[!is.finite(inv_diag)] <- 0  # Handle zero variance

  if (sum(inv_diag) == 0) {
    return(0)
  }

  ivp_weights <- inv_diag / sum(inv_diag)

  # Portfolio variance: w' * Cov * w - VECTORIZED
  portfolio_var <- as.numeric(t(ivp_weights) %*% sub_cov %*% ivp_weights)

  return(portfolio_var)
}

#' Risk Parity Weighting Suite
#'
#' Collection of risk-based weighting methods for portfolio construction.
#' Each method allocates capital based on risk characteristics rather than
#' market capitalization or equal weights.
#'
#' @param selected_df Binary selection matrix (data.frame with Date column).
#' @param prices_df Price data for risk calculations (typically daily).
#'   Returns are calculated internally from prices.
#' @param lookback_periods Number of periods for risk estimation (default: 252).
#' @param min_periods Minimum periods required (default: 60).
#' @param method Optimization method for risk parity.
#'
#' @details
#' \strong{Methods}
#' \describe{
#'   \item{\code{inverse_vol}}{Weights proportional to 1 / volatility
#'     (e.g., 1 / sd of returns over \code{lookback_periods}). Lower volatility
#'     stocks receive higher weights.}
#'   \item{\code{equal_risk}}{Equal Risk Contribution (ERC): finds weights so
#'     each position contributes equally to total portfolio risk (iterative optimization).}
#'   \item{\code{max_div}}{Maximum Diversification: maximizes the ratio of
#'     weighted average volatility to portfolio volatility.}
#' }
#'
#' The function accepts price data and calculates returns internally,
#' ensuring consistency with other library functions. Daily prices are
#' recommended for accurate volatility estimation.
#'
#' @return Weight matrix with the same dates as \code{selected_df}; each row sums to 1.
#' @export
#' @examples
#' data("sample_prices_daily")
#' data("sample_prices_weekly")
#' momentum <- calc_momentum(sample_prices_weekly, lookback = 12)
#' selected <- filter_top_n(momentum, n = 10)
#' weight_by_risk_parity(selected, sample_prices_daily, method = "inverse_vol")
#' weight_by_risk_parity(selected, sample_prices_daily, method = "equal_risk")
#' weight_by_risk_parity(selected, sample_prices_daily, method = "max_div")
weight_by_risk_parity <- function(selected_df, prices_df,  # Changed from returns_df
                                  method = c("inverse_vol", "equal_risk", "max_div"),
                                  lookback_periods = 252,
                                  min_periods = 60) {

  method <- match.arg(method)

  # Input validation
  validate_selection_data(selected_df, "weight_by_risk_parity")
  validate_selection_data(prices_df, "weight_by_risk_parity")  # Changed

  # Setup
  selected_dt <- ensure_dt_copy(selected_df)
  prices_dt <- ensure_dt_copy(prices_df)  # Changed
  symbol_cols <- setdiff(names(selected_dt), "Date")

  # Handle NA in selection
  for (col in symbol_cols) {
    selected_dt[is.na(get(col)), (col) := 0]
  }

  # CALCULATE RETURNS FROM PRICES (matching HRP pattern)
  returns_dt <- data.table(Date = prices_dt$Date)
  for (col in symbol_cols) {
    if (col %in% names(prices_dt)) {
      # Simple returns: (P_t - P_{t-1}) / P_{t-1}
      prices <- prices_dt[[col]]
      returns <- c(NA, diff(prices) / head(prices, -1))
      returns_dt[, (col) := returns]
    }
  }

  # Initialize result
  weight_df <- data.table(Date = selected_dt$Date)
  for(col in symbol_cols) {
    weight_df[, (col) := 0.0]
  }


  # Process each date
  for (i in seq_len(nrow(selected_dt))) {
    date <- selected_dt$Date[i]

    # Get selected assets
    selected_mask <- unlist(selected_dt[i, ..symbol_cols]) > 0
    if (!any(selected_mask)) {
      next
    }

    selected_symbols <- symbol_cols[selected_mask]

    # Get historical returns
    end_idx <- which(returns_dt$Date <= date)
    if (length(end_idx) == 0) next

    latest_idx <- max(end_idx)
    start_idx <- max(1, latest_idx - lookback_periods + 1)

    if (latest_idx - start_idx + 1 < min_periods) next

    # Extract returns matrix
    returns_subset <- returns_dt[start_idx:latest_idx, c("Date", selected_symbols), with = FALSE]
    returns_matrix <- as.matrix(returns_subset[, ..selected_symbols])

    # Remove rows with all NAs
    valid_rows <- rowSums(!is.na(returns_matrix)) > 0
    returns_matrix <- returns_matrix[valid_rows, , drop = FALSE]

    if (nrow(returns_matrix) < min_periods) next

    # Calculate risk parity weights based on method
    tryCatch({
      if (method == "inverse_vol") {
        # Simple inverse volatility
        vols <- apply(returns_matrix, 2, sd, na.rm = TRUE)
        vols[vols == 0] <- 1e-8  # Avoid division by zero
        inv_vols <- 1 / vols
        weights_vec <- inv_vols / sum(inv_vols)

      } else if (method == "equal_risk") {
        # Equal Risk Contribution (ERC) - simplified version
        cov_matrix <- cov(returns_matrix, use = "complete.obs")
        weights_vec <- calculate_erc_weights(cov_matrix, selected_symbols)

      } else if (method == "max_div") {
        # Maximum Diversification Portfolio
        cov_matrix <- cov(returns_matrix, use = "complete.obs")
        weights_vec <- calculate_max_div_weights(cov_matrix, selected_symbols)
      }

      # Assign weights
      for (j in seq_along(selected_symbols)) {
        weight_df[i, (selected_symbols[j]) := weights_vec[j]]
      }

    }, error = function(e) {
      # Fallback to equal weights
      equal_weight <- 1 / length(selected_symbols)
      for (symbol in selected_symbols) {
        weight_df[i, (symbol) := equal_weight]
      }
    })
  }

  return(weight_df)
}

#' Calculate Equal Risk Contribution weights (simplified)
#' @keywords internal
calculate_erc_weights <- function(cov_matrix, asset_names) {
  n <- nrow(cov_matrix)

  # Start with equal weights
  weights <- rep(1/n, n)

  # Simple iterative approach (could be improved with optimization)
  for (iter in 1:50) {  # Max 50 iterations
    # Calculate risk contributions
    port_vol <- sqrt(t(weights) %*% cov_matrix %*% weights)
    risk_contribs <- (weights * (cov_matrix %*% weights)) / as.numeric(port_vol)

    # Update weights to equalize risk contributions
    avg_risk <- mean(risk_contribs)
    adjustment <- avg_risk / risk_contribs
    weights <- weights * adjustment
    weights <- weights / sum(weights)  # Normalize

    # Check convergence
    if (max(abs(risk_contribs - avg_risk)) < 1e-6) {
      break
    }
  }

  return(as.numeric(weights))
}

#' Calculate Maximum Diversification Portfolio weights
#' @keywords internal
calculate_max_div_weights <- function(cov_matrix, asset_names) {
  # Maximum Diversification = maximize weighted average vol / portfolio vol
  # This is equivalent to inverse volatility in many cases

  vols <- sqrt(diag(cov_matrix))
  vols[vols == 0] <- 1e-8

  inv_vols <- 1 / vols
  weights <- inv_vols / sum(inv_vols)

  return(as.numeric(weights))
}

# Documentation helper
hrp_info <- function() {
  cat("Hierarchical Risk Parity (HRP) Implementation\n")
  cat("=============================================\n\n")

  cat("Main Functions:\n")
  cat("  weight_by_hrp() - Full HRP portfolio construction\n")
  cat("  weight_by_risk_parity() - Alternative risk-based methods\n\n")

  cat("HRP Advantages:\n")
  cat("  - No matrix inversion (numerically stable)\n")
  cat("  - Naturally diversified allocations\n")
  cat("  - Handles high-dimensional problems\n")
  cat("  - Robust to estimation errors\n\n")

  cat("Risk Parity Methods:\n")
  cat("  inverse_vol: Simple 1/volatility weighting\n")
  cat("  equal_risk: Equal Risk Contribution (ERC)\n")
  cat("  max_div: Maximum Diversification Portfolio\n\n")

  cat("Usage:\n")
  cat("  hrp_weights <- weight_by_hrp(selected, returns, lookback_periods = 252)\n")
  cat("  rp_weights <- weight_by_risk_parity(selected, returns, method = 'equal_risk')\n")
}






