# utils.R
# Utility functions for safe operations in momentum strategies library
# These functions prevent common errors with NA values, zero division, etc.

###############################################################################
# SAFE ARITHMETIC OPERATIONS
###############################################################################

#' Safe Division with NA and Zero Handling
#'
#' @description
#' Performs division with automatic handling of NA values, zeros, and infinity.
#' Returns 0 for division by zero and NA cases.
#'
#' @param numerator Numeric vector
#' @param denominator Numeric vector
#'
#' @return Numeric vector with safe division results
#' @export
#' @examples
#' safe_divide(c(10, 0, NA, 5), c(2, 0, 5, NA))  # Returns c(5, 0, 0, 0)
safe_divide <- function(numerator, denominator) {
  result <- numerator / denominator
  result[is.na(result) | is.infinite(result)] <- 0
  return(result)
}
###############################################################################
# SAFE LOGICAL OPERATIONS
###############################################################################

#' Safe ANY Operation with NA Handling
#'
#' @description
#' Performs any() operation that returns FALSE when all values are NA
#' instead of NA.
#'
#' @param x Logical vector
#' @param ... Additional arguments passed to any()
#'
#' @return Logical value (never NA)
#' @keywords internal
safe_any <- function(x, ...) {
  if (all(is.na(x))) return(FALSE)
  return(any(x, na.rm = TRUE, ...))
}

# Safe logical AND that handles NAs
safe_and <- function(x, y) {
  result <- rep(FALSE, length(x))
  valid_idx <- !is.na(x) & !is.na(y)
  result[valid_idx] <- x[valid_idx] & y[valid_idx]
  return(result)
}

###############################################################################
# DATA TYPE OPERATIONS
###############################################################################

#' Ensure Data.Table Without Mutation
#'
#' @description
#' Converts input to data.table if needed, always returning a copy
#' to prevent accidental data mutation. Core safety function used
#' throughout the library.
#'
#' @param data Data.frame or data.table
#'
#' @return Copy of data as data.table
#' @export
#' @examples
#' data("sample_prices_weekly")
#' dt <- ensure_dt_copy(sample_prices_weekly)  # Safe to modify dt
ensure_dt_copy <- function(data) {
  if (!is.data.table(data)) {
    return(as.data.table(data))
  } else {
    return(copy(data))
  }
}

###############################################################################
# VALIDATION HELPERS
###############################################################################

# Check if all values are valid (not NA or Inf)
all_valid <- function(x) {
  all(!is.na(x) & !is.infinite(x))
}

# Count valid (non-NA) values
count_valid <- function(x) {
  sum(!is.na(x))
}

###############################################################################
# ADDITIONAL SAFETY FUNCTIONS
###############################################################################

# Safe sum that returns 0 if all values are NA
safe_sum <- function(x, ...) {
  if (all(is.na(x))) return(0)
  return(sum(x, na.rm = TRUE, ...))
}

# Safe mean that returns NA if all values are NA
safe_mean <- function(x, ...) {
  if (all(is.na(x))) return(NA_real_)
  return(mean(x, na.rm = TRUE, ...))
}

# Safe max that returns -Inf if all values are NA
safe_max <- function(x, ...) {
  if (all(is.na(x))) return(-Inf)
  return(max(x, na.rm = TRUE, ...))
}

# Safe min that returns Inf if all values are NA
safe_min <- function(x, ...) {
  if (all(is.na(x))) return(Inf)
  return(min(x, na.rm = TRUE, ...))
}

###############################################################################
# STRATEGY HELPER FUNCTIONS
###############################################################################

#' Convert Conditions to Selection Format
#'
#' @description
#' Converts condition matrices or data frames to standard selection format
#' with Date column and binary values. Handles NA by converting to 0.
#'
#' @param condition_matrix Matrix or data frame with conditions
#' @param date_column Optional Date vector if not in input
#'
#' @return Data.table in selection format (Date + binary columns)
#' @export
#' @examples
#' data("sample_prices_weekly")
#' ma20 <- calc_moving_average(sample_prices_weekly, 20)
#' above_ma <- filter_above(calc_distance(sample_prices_weekly, ma20), 0)
#' selection <- as_selection(above_ma, sample_prices_weekly$Date)
as_selection <- function(condition_matrix, date_column = NULL) {
  # Handle different input types
  if (is.data.frame(condition_matrix) && "Date" %in% names(condition_matrix)) {
    # Already has Date column
    result <- as.data.table(condition_matrix)
  } else if (is.matrix(condition_matrix) || is.data.frame(condition_matrix)) {
    # Convert to data.table and add Date
    result <- as.data.table(condition_matrix)
    if (!is.null(date_column)) {
      result$Date <- date_column
    } else {
      stop("Need to provide date_column for matrix input")
    }
  } else {
    stop("Input must be a matrix or data.frame")
  }

  # Convert logical to numeric and handle NA
  symbol_cols <- setdiff(names(result), "Date")
  for (col in symbol_cols) {
    # Convert logical to numeric
    if (is.logical(result[[col]])) {
      result[, (col) := as.numeric(get(col))]
    }
    # Convert NA to 0
    result[is.na(get(col)), (col) := 0]
  }

  return(result)
}

# Run complete strategy with natural condition syntax
run_strategy <- function(prices, selection_condition,
                         weighting_method = "equal",
                         signal_for_weighting = NULL,
                         initial_capital = 100000,
                         name = "Strategy",
                         skip_warmup = TRUE,
                         verbose = FALSE) {

  # Convert natural condition to selection format
  selected <- as_selection(selection_condition, prices$Date)

  # Optional: Skip warmup period
  if (skip_warmup) {
    symbol_cols <- setdiff(names(selected), "Date")
    # Find first row where at least one stock is selected
    first_valid <- which(rowSums(selected[, ..symbol_cols]) > 0)[1]

    if (!is.na(first_valid) && first_valid > 1 && verbose) {
      message(sprintf("Skipping %d warmup periods", first_valid - 1))
    }
  }

  # Apply weighting
  if (weighting_method == "equal") {
    weights <- weight_equally(selected)
  } else if (weighting_method == "signal" && !is.null(signal_for_weighting)) {
    signal_selection <- as_selection(signal_for_weighting, prices$Date)
    weights <- weight_by_signal(selected, signal_selection)
  } else if (weighting_method == "rank" && !is.null(signal_for_weighting)) {
    signal_selection <- as_selection(signal_for_weighting, prices$Date)
    weights <- weight_by_rank(selected, signal_selection)
  } else if (weighting_method == "volatility") {
    # For volatility weighting, we need additional parameters
    # This is a simplified version - you might want to expand this
    weights <- weight_by_volatility(selected, prices, low_vol_preference = TRUE)
  } else {
    stop("Invalid weighting method. Choose: 'equal', 'signal', 'rank', or 'volatility'")
  }

  # Run backtest
  run_backtest(prices, weights, initial_capital, name, verbose = verbose)
}
###############################################################################
# REGIME FILTER FUNCTIONS
###############################################################################

# Internal helper - efficient regime application
apply_regime_broadcast <- function(selection_matrix, regime_condition) {
  result <- ensure_dt_copy(selection_matrix)
  symbol_cols <- setdiff(names(result), "Date")

  regime_binary <- as.numeric(!is.na(regime_condition) & regime_condition)

  for (col in symbol_cols) {
    result[, (col) := get(col) * regime_binary]
  }

  return(result)
}

# Main user function - filter selections based on market regime
filter_when <- function(selection_matrix, condition, label = NULL) {
  result <- apply_regime_broadcast(selection_matrix, condition)

  if (!is.null(label) && interactive()) {
    active_pct <- sum(condition, na.rm = TRUE) / sum(!is.na(condition)) * 100
    message(sprintf("%s filter: active %.1f%% of time", label, active_pct))
  }

  return(result)
}


###############################################################################
# DOCUMENTATION
###############################################################################

# Print utility functions information
utils_info <- function() {
  cat("Momentum Strategies Library - Utility Functions\n")
  cat("==============================================\n\n")

  cat("Safe Arithmetic:\n")
  cat("  safe_divide(x, y) - Division with zero/NA protection\n\n")

  cat("Safe Logical:\n")
  cat("  safe_any(x) - any() that handles all NA\n")
  cat("  safe_and(x, y) - Logical AND with NA handling\n\n")

  cat("Data Operations:\n")
  cat("  ensure_dt_copy(data) - Ensure data.table without mutation\n\n")

  cat("Validation:\n")
  cat("  all_valid(x) - Check for NA/Inf values\n")
  cat("  count_valid(x) - Count non-NA values\n\n")

  cat("Additional:\n")
  cat("  safe_sum(), safe_mean(), safe_max(), safe_min()\n\n")

  cat("Strategy Helpers:\n")
  cat("  as_selection(condition, dates) - Convert conditions to selection format\n")
  cat("  run_strategy(...) - Complete strategy runner with natural syntax\n")
}





# Helper for NULL coalescing (add to utils.R if not exists)
`%||%` <- function(x, y) if (is.null(x)) y else x



#' Align Data to Strategy Timeframe
#'
#' @description
#' Aligns higher-frequency data to match strategy timeframe.
#'
#' @param high_freq_data Data frame to align
#' @param low_freq_dates Date vector from strategy
#' @param method Alignment method: "forward_fill", "nearest", or "interpolate"
#'
#' @return Aligned data frame
#' @export
#' @examples
#' data("sample_prices_weekly")
#' data("sample_prices_daily")
#' momentum <- calc_momentum(sample_prices_weekly, lookback = 12)
#' selected <- filter_top_n(momentum, 10)
#' # Create a stability signal from daily data
#' daily_vol <- calc_rolling_volatility(sample_prices_daily, lookback = 20)
#' stability_signal <- align_to_timeframe(daily_vol, sample_prices_weekly$Date)
#' weights <- weight_by_signal(selected, stability_signal)
align_to_timeframe <- function(high_freq_data, low_freq_dates, method = c("forward_fill", "nearest", "interpolate")) {
  # Align high-frequency data to low-frequency dates
  # OPTIMIZED VERSION using data.table rolling joins for speed
  #
  # This function solves the common problem of aligning indicators calculated
  # on one timeframe (e.g., daily) to a trading timeframe (e.g., weekly).
  # For example, align daily volatility calculations to weekly trading dates.
  #
  # Args:
  #   high_freq_data: data.frame/data.table with Date column and data columns
  #                   Example: daily prices or daily indicators
  #   low_freq_dates: Vector of dates to align to
  #                   Example: weekly trading dates
  #   method: Alignment method
  #     - "forward_fill": Use most recent available value (default)
  #                       Best for: prices, cumulative indicators
  #     - "nearest": Use value from nearest date (before or after)
  #                  Best for: volatility, correlations
  #     - "interpolate": Linear interpolation between values
  #                      Best for: smooth indicators like moving averages
  #
  # Returns:
  #   data.table with Date column matching low_freq_dates and all data
  #   columns from high_freq_data aligned according to method
  #
  # Examples:
  #   # Align daily volatility to weekly dates
  #   daily_vol <- calc_rolling_volatility(daily_prices, 20)
  #   weekly_vol <- align_to_timeframe(daily_vol, weekly_dates, "forward_fill")
  #
  #   # Align daily correlation for weekly strategy
  #   daily_corr <- calc_rolling_correlation(daily_prices, "SPY", 60)
  #   weekly_corr <- align_to_timeframe(daily_corr, weekly_dates, "forward_fill")

  method <- match.arg(method)

  # Ensure data.table for fast operations
  dt <- ensure_dt_copy(high_freq_data)
  setkey(dt, Date)  # Set key for ultra-fast date joins

  # Get column names (excluding Date)
  data_cols <- setdiff(names(dt), "Date")

  # Create target dates table
  target_dt <- data.table(Date = low_freq_dates)
  setkey(target_dt, Date)

  if (method == "forward_fill") {
    # Use data.table's rolling join with roll = TRUE
    # This finds the most recent date that is <= target date
    # MUCH faster than looping through dates
    result <- dt[target_dt, roll = TRUE]

  } else if (method == "nearest") {
    # Use data.table's rolling join with roll = "nearest"
    # This finds the closest date (either before or after)
    result <- dt[target_dt, roll = "nearest"]

  } else {  # interpolate
    # For interpolation, use R's approx() function
    # Already efficient for this use case
    result <- copy(target_dt)

    # Interpolate all columns
    for (col in data_cols) {
      interp_values <- approx(
        x = as.numeric(dt$Date),
        y = dt[[col]],
        xout = as.numeric(low_freq_dates),
        method = "linear",
        rule = 2  # Use nearest value for extrapolation beyond data range
      )$y

      result[, (col) := interp_values]
    }
  }

  # Ensure the result has exactly the requested dates
  # (important because rolling joins might have different dates)
  if (method != "interpolate") {
    result[, Date := low_freq_dates]
  }

  return(result)
}






#' Invert Signal Values for Preference Reversal
#'
#' @description
#' Transforms signal values using (1 - value) to reverse preference direction.
#' Useful when high values indicate something to avoid. For example, inverting
#' volatility makes low-vol stocks appear as high signals.
#'
#' @param signal_df Data frame with Date column and signal columns
#'
#' @return Data frame with inverted signal values
#' @export
#' @examples
#' data("sample_prices_weekly")
#' # Prefer low volatility stocks
#' volatility <- calc_rolling_volatility(sample_prices_weekly, 20)
#' stability_signal <- invert_signal(volatility)
#' # Select top 10 momentum stocks first
#' momentum <- calc_momentum(sample_prices_weekly, 12)
#' selected <- filter_top_n(momentum, 10)
#' # Weight by inverted volatility (low vol = high weight)
#' weights <- weight_by_signal(selected, stability_signal)
invert_signal <- function(signal_df) {
  # Input validation
  if (!is.data.frame(signal_df)) {
    stop("invert_signal: signal_df must be a data.frame or data.table")
  }

  if (!"Date" %in% names(signal_df)) {
    stop("invert_signal: signal_df must have a 'Date' column")
  }

  if (ncol(signal_df) < 2) {
    stop("invert_signal: signal_df must have at least one signal column besides Date")
  }

  # Create copy to avoid modifying input
  dt <- ensure_dt_copy(signal_df)

  # Get signal columns (all except Date)
  signal_cols <- setdiff(names(dt), "Date")

  # Apply inversion: new_value = 1 - old_value
  # This works for any numeric values, though normalized [0,1] is ideal
  dt[, (signal_cols) := lapply(.SD, function(x) 1 - x), .SDcols = signal_cols]

  return(dt)
}

#' Internal helper: simple coalesce
#' @keywords internal
#' @noRd
.pt_coalesce <- function(x, default) {
  if (is.null(x) || (length(x) == 1 && is.na(x))) default else x
}

#' Internal helper: coalesce limited to NULL or NA
#' @keywords internal
#' @noRd
.pt_null_coalesce <- function(x, default) .pt_coalesce(x, default)

#' Internal helper: detect sampling frequency (per year)
#' @keywords internal
#' @noRd
.pt_detect_frequency <- function(dates) {
  if (length(dates) < 2) return(252)
  d <- as.Date(dates)
  diffs <- as.numeric(diff(sort(d)))
  med   <- suppressWarnings(median(diffs, na.rm = TRUE))
  if (!is.finite(med)) return(252)

  if (med <= 1.5) return(252)                    # daily
  if (med >= 6   && med <= 8)  return(52)        # weekly
  if (med >= 28  && med <= 32) return(12)        # monthly
  if (med >= 85  && med <= 95) return(4)         # quarterly

  span <- as.numeric(max(d) - min(d))
  if (is.finite(span) && span > 0) {
    ppy <- round(length(d) * 365.25 / span)
    if (ppy > 100) return(252)
    if (ppy >  40) return(52)
    if (ppy >   8) return(12)
    return(4)
  }
  252
}

#' Internal helper: reject cadence/timeframe knobs in grids
#' @keywords internal
#' @noRd
.pt_check_grid_params <- function(grid, context = "optimization") {
  # Block structural fields + cadence/timeframe synonyms
  forbidden <- c(
    # structural / reserved
    "date","prices","returns","weights",
    # cadence / timeframe
    "rebalance","rebalance_period","rebalance_freq","rebalance_by",
    "rebalance_every","cadence","timeframe","ta_timeframe",
    "trade_freq","trade_period","frequency","rebal_weeks"
  )
  nm  <- tolower(names(grid))
  bad <- intersect(nm, forbidden)
  if (length(bad)) {
    stop(sprintf(
      "Grid contains forbidden parameter names for %s: %s",
      context, paste(bad, collapse = ", ")
    ), call. = FALSE)
  }
}

#' Internal helper: robust weights validator and aligner
#' @keywords internal
#' @noRd
.pt_validate_weights <- function(prices, weights, context = "weights", allow_sparse = TRUE) {
  # Work internally as base data.frames
  P <- as.data.frame(prices,  check.names = FALSE, stringsAsFactors = FALSE)
  W <- as.data.frame(weights, check.names = FALSE, stringsAsFactors = FALSE)

  if (!"Date" %in% names(P)) stop("prices must have a 'Date' column", call. = FALSE)
  if (!"Date" %in% names(W)) stop(sprintf("[%s] weights must have a 'Date' column", context), call. = FALSE)

  # 1) Normalize Date
  to_date <- function(x) {
    if (inherits(x, "Date"))   return(x)
    if (inherits(x, "POSIXt")) return(as.Date(x))
    if (is.numeric(x))         return(as.Date(x, origin = "1970-01-01"))
    as.Date(x)
  }
  P$Date <- to_date(P$Date)
  W$Date <- to_date(W$Date)

  # Collapse duplicate dates in weights (keep last)
  if (anyDuplicated(W$Date)) {
    W <- W[order(W$Date), , drop = FALSE]
    W <- W[!duplicated(W$Date, fromLast = TRUE), , drop = FALSE]
    warning(sprintf("[%s] Duplicate dates in weights; keeping last per date.", context), call. = FALSE)
  }

  # 2) Columns: drop junk, align to prices
  asset_p <- setdiff(names(P), "Date")
  asset_w <- setdiff(names(W), "Date")

  extra <- setdiff(asset_w, asset_p)
  if (length(extra)) {
    warning(sprintf("[%s] Dropping extra non-asset columns: %s",
                    context, paste(extra, collapse = ", ")), call. = FALSE)
    asset_w <- intersect(asset_w, asset_p)
    W <- W[, c("Date", asset_w), drop = FALSE]
  }
  if (!length(asset_w)) {
    warning(sprintf("[%s] No matching asset columns; returning all-zero weights.", context), call. = FALSE)
    out_df <- data.frame(Date = P$Date, check.names = FALSE)
    for (nm in asset_p) out_df[[nm]] <- 0
    # Return as data.table if available (so run_backtest's with=FALSE works)
    if (requireNamespace("data.table", quietly = TRUE)) {
      return(data.table::as.data.table(out_df))
    } else {
      return(out_df)
    }
  }

  common_assets <- intersect(asset_p, asset_w)
  W <- W[, c("Date", common_assets), drop = FALSE]

  # 3) Date alignment + diagnostics
  idx <- match(P$Date, W$Date)
  present <- !is.na(idx)
  match_rate <- mean(present)

  if (!all(present) && !allow_sparse) {
    stop(sprintf("[%s] weights missing %d/%d required dates (allow_sparse=FALSE).",
                 context, sum(!present), length(idx)), call. = FALSE)
  }
  if (match_rate < 0.99) {
    warning(sprintf("[%s] Only %.1f%% of dates matched; non-matches set to 0. Check Date class/timezone.",
                    context, 100 * match_rate), call. = FALSE)
  }

  # Build aligned output (zeros by default)
  out_df <- data.frame(Date = P$Date, check.names = FALSE)
  for (nm in asset_p) out_df[[nm]] <- 0
  if (any(present)) for (nm in common_assets) out_df[present, nm] <- W[idx[present], nm]

  # 4) Coerce numeric, NA->0, trim tiny noise
  for (nm in asset_p) {
    out_df[[nm]] <- suppressWarnings(as.numeric(out_df[[nm]]))
    if (anyNA(out_df[[nm]])) {
      na_n <- sum(is.na(out_df[[nm]]))
      warning(sprintf("[%s] NAs in weights for %s (%d); replacing with 0.", context, nm, na_n), call. = FALSE)
      out_df[[nm]][is.na(out_df[[nm]])] <- 0
    }
    out_df[[nm]][abs(out_df[[nm]]) < 1e-14] <- 0
  }

  # 5) Degeneracy hint
  tradable <- sum(rowSums(abs(out_df[, asset_p, drop = FALSE])) > 0)
  if (tradable < 2L) {
    warning(sprintf("[%s] Weights are effectively all-zero on almost all dates after alignment (%d/%d tradable rows). Check lookback/selection/date alignment.",
                    context, tradable, nrow(out_df)), call. = FALSE)
  }

  # Return as data.table if available (keeps backtester happy)
  if (requireNamespace("data.table", quietly = TRUE)) {
    out_df <- data.table::as.data.table(out_df)
  }
  out_df
}

#' Volatility targeting (row-wise) with optional down-only cap
#'
#' Scales each row of a wide weight table (`Date + symbols`) so the estimated
#' annualised portfolio volatility matches a target. Volatility is estimated
#' from a rolling covariance of simple asset returns computed from `prices`.
#'
#' Weights decided at t-1 apply to returns over t.
#'
#' @param weights data.frame/data.table with columns: `Date`, then one column per symbol.
#' @param prices data.frame/data.table of adjusted prices at the same cadence as `weights`:
#'   first column `Date`, remaining columns one per symbol (matching weight symbols).
#' @param lookback Integer window length (in periods) for the rolling covariance. Default 26.
#' @param target_annual Annualised volatility target (e.g., 0.12 for 12%). Must be > 0.
#' @param periods_per_year Number of periods per year used for annualisation (e.g., 52 for weekly).
#' @param cap If TRUE (default), scale down only: exposure is reduced when the
#'   estimated vol exceeds the target, and untouched otherwise. In this down-only mode,
#'   the function adds/overwrites a `CASH` column equal to
#'   `1 - rowSums(pmax(scaled_weights, 0))` so that symbols + CASH is approximately 1.
#'   If FALSE, the scaler may be > 1 (leverage allowed) and no `CASH` is added.
#'
#' @return A `data.table` with the same `Date` and symbol columns as `weights`
#'   (plus `CASH` when `cap = TRUE`).
#'
#' @details
#' The covariance at row i is computed from the last `lookback` rows of simple
#' returns up to i (inclusive), estimated on the intersection of symbols present
#' in both `weights` and `prices`. The row scaler is
#' `s_i = min(1, target_vol / est_vol_i)` when `cap = TRUE`, and
#' `s_i = target_vol / est_vol_i` when `cap = FALSE`, with safeguards for zero or
#' non-finite variances.
#'
#' @examples
#' \donttest{
#'   data(sample_prices_weekly)
#'   mom12 <- PortfolioTesteR::calc_momentum(sample_prices_weekly, 12)
#'   sel10 <- PortfolioTesteR::filter_top_n(mom12, 10)
#'   w_eq  <- PortfolioTesteR::weight_equally(sel10)
#'
#'   w_vt  <- vol_target(w_eq, sample_prices_weekly, lookback = 26,
#'                       target_annual = 0.12, periods_per_year = 52, cap = TRUE)
#'   head(w_vt)
#' }
#'
#' @export
vol_target <- function(weights,
                       prices,
                       lookback = 26L,
                       target_annual = 0.12,
                       periods_per_year = 52L,
                       cap = TRUE) {
  CASH <- NULL  # for data.table NSE, silences R CMD check
  stopifnot(is.data.frame(weights), is.data.frame(prices))
  stopifnot("Date" %in% names(weights), "Date" %in% names(prices))
  stopifnot(is.numeric(lookback), lookback >= 2L)
  stopifnot(is.numeric(target_annual), is.finite(target_annual), target_annual > 0)
  stopifnot(is.numeric(periods_per_year), periods_per_year > 0)

  W <- data.table::as.data.table(data.table::copy(weights))
  P <- data.table::as.data.table(data.table::copy(prices))
  data.table::setkey(W, Date); data.table::setkey(P, Date)

  # intersect symbols that exist in both
  syms_w <- setdiff(names(W), "Date")
  syms_p <- setdiff(names(P), "Date")
  syms   <- intersect(syms_w, syms_p)
  if (!length(syms)) stop("vol_target(): no overlapping symbols between weights and prices.")

  # simple returns (Date + syms); returns over t use prices[t-1], prices[t]
  R <- PortfolioTesteR::panel_returns_simple(P[, c("Date", syms), with = FALSE])

  # Align shapes for decision lag: weights at t-1 apply to returns at t
  # We'll compute a scaler for each row of W using returns up to that decision time.
  W_sym <- as.matrix(W[, ..syms])
  nW    <- nrow(W)
  nR    <- nrow(R)

  # Pre-compute simple returns matrix (drop Date)
  R_mat <- as.matrix(R[, ..syms])

  # Rolling variance of portfolio each row i (using last 'lookback' rows of returns up to i)
  # Note: R has one fewer row than prices; we clamp indices safely.
  scaler <- rep(1, nW)
  for (i in seq_len(nW)) {
    # index in R that corresponds to decision time i is min(i, nR)
    ri <- min(i, nR)
    if (ri >= lookback) {
      idx <- (ri - lookback + 1L):ri
      S   <- try(stats::cov(R_mat[idx, , drop = FALSE], use = "pairwise.complete.obs"),
                 silent = TRUE)
      if (inherits(S, "try-error") || anyNA(S)) {
        s_i <- 1
      } else {
        w_i <- W_sym[i, ]
        # guard: missing weights treated as 0
        w_i[!is.finite(w_i)] <- 0
        v   <- as.numeric(drop(t(w_i) %*% S %*% w_i))
        sig <- sqrt(max(v, 0)) * sqrt(periods_per_year)
        if (!is.finite(sig) || sig == 0) {
          s_i <- 1
        } else {
          s_i <- target_annual / sig
          if (isTRUE(cap)) s_i <- min(1, s_i)
        }
      }
      scaler[i] <- s_i
    } else {
      scaler[i] <- 1
    }
  }

  # Apply scaler row-wise to symbol weights
  W_scaled <- data.table::copy(W)
  W_scaled[, (syms) := as.data.table(W_sym * scaler)]

  # When down-only: add/overwrite CASH so (symbols + CASH) approximately 1 (long side budget)
  if (isTRUE(cap)) {
    pos_sum <- rowSums(pmax(as.matrix(W_scaled[, ..syms]), 0), na.rm = TRUE)
    cash    <- 1 - pos_sum
    cash[cash < 0 | !is.finite(cash)] <- 0
    W_scaled[, CASH := cash]
  }

  W_scaled[]
}

#' Cap per-symbol and per-group exposures, with optional renormalization
#'
#' Row-wise caps on a wide weights table (`Date` + symbols):
#' \enumerate{
#'   \item per-symbol cap (`max_per_symbol`)
#'   \item optional per-group cap (`max_per_group`) using `group_map`
#' }
#'
#' Negatives are clamped to 0 unless `allow_short = TRUE`.
#'
#' Renormalization options:
#' \itemize{
#'   \item \code{renormalize = "none"}: leave gross (sum of positive weights) as is.
#'   \item \code{renormalize = "down"}: if gross > 1, scale down so gross == 1.
#'   \item \code{renormalize = "both"}: if 0 < gross != 1, scale so gross == 1.
#' }
#'
#' The argument \code{renorm} (logical) is a deprecated alias;
#' \code{TRUE} behaves like \code{renormalize = "both"}.
#'
#' The \code{caps} list form is supported:
#' \describe{
#'   \item{max_per_symbol}{Maximum absolute weight per symbol (0-1).}
#'   \item{max_per_group}{Maximum gross weight per group (0-1).}
#'   \item{group_map}{Named character vector mapping \code{symbol -> group}.}
#' }
#'
#' @param weights data.frame/data.table: columns `Date`, then symbols.
#' @param max_per_symbol numeric scalar or named vector (absolute cap per symbol).
#' @param group_map optional named character vector or data.frame(`symbol`,`group`).
#' @param max_per_group optional numeric scalar (per-group gross cap).
#' @param allow_short logical; if `FALSE` clamp negatives to 0.
#' @param renormalize character: one of `"none"`, `"down"`, `"both"`. Default `"none"`.
#' @param renormalize_down logical; deprecated alias for down-only (kept for compatibility).
#' @param renorm logical; deprecated alias; if `TRUE` acts like `renormalize="both"`.
#' @param cash_col optional character; if provided, set to `1 - sum(pmax(w,0))`.
#' @param caps optional list; alternative to split args (see details above).
#' @return data.table of capped (and optionally renormalized) weights.
#' @export
cap_exposure <- function(weights,
                         max_per_symbol    = NULL,
                         group_map         = NULL,
                         max_per_group     = NULL,
                         allow_short       = FALSE,
                         renormalize       = c("none","down","both"),
                         renormalize_down  = FALSE,
                         renorm            = NULL,
                         cash_col          = NULL,
                         caps              = NULL)
{
  stopifnot(is.data.frame(weights), "Date" %in% names(weights))
  DT <- data.table::as.data.table(data.table::copy(weights))
  data.table::setkey(DT, Date)
  syms <- setdiff(names(DT), "Date")
  if (!length(syms)) stop("cap_exposure(): no symbol columns found.")

  # ----- accept list-form 'caps' -----
  if (!is.null(caps)) {
    if (is.null(max_per_symbol) && !is.null(caps$max_per_symbol)) max_per_symbol <- caps$max_per_symbol
    if (is.null(max_per_group)  && !is.null(caps$max_per_group))  max_per_group  <- caps$max_per_group
    if (is.null(group_map)      && !is.null(caps$group_map))      group_map      <- caps$group_map
  }

  # ----- per-symbol cap vector -----
  cap_sym <- rep_len(as.numeric(if (is.null(max_per_symbol)) Inf else max_per_symbol), length(syms))
  names(cap_sym) <- syms
  if (!is.null(max_per_symbol) && length(max_per_symbol) > 1L) {
    if (!is.null(names(max_per_symbol))) {
      cap_sym[names(max_per_symbol)] <- as.numeric(max_per_symbol)
    } else {
      cap_sym <- as.numeric(max_per_symbol)
      if (length(cap_sym) != length(syms))
        stop("cap_exposure(): max_per_symbol length must be 1 or match #symbols.")
      names(cap_sym) <- syms
    }
  }

  # ----- group map -----
  grp <- NULL
  if (!is.null(group_map) && !is.null(max_per_group)) {
    if (is.data.frame(group_map)) {
      stopifnot(all(c("symbol","group") %in% names(group_map)))
      grp <- setNames(as.character(group_map$group), group_map$symbol)
    } else {
      if (is.null(names(group_map)))
        stop("cap_exposure(): group_map must be named (names = symbols).")
      grp <- as.character(group_map)
      names(grp) <- names(group_map)
    }
    grp <- grp[intersect(names(grp), syms)]
    if (!length(grp)) grp <- NULL
  }

  # ----- pick renormalize mode (honour legacy flags) -----
  renormalize <- match.arg(renormalize)
  if (isTRUE(renormalize_down) && identical(renormalize, "none")) renormalize <- "down"
  if (!is.null(renorm)) {
    if (isTRUE(renorm)) renormalize <- "both"
  }

  # helpers
  clamp_symbol <- function(w_row) {
    w <- as.numeric(w_row)
    if (!allow_short) w[w < 0] <- 0
    caps <- cap_sym
    if (allow_short) {
      w <- pmax(pmin(w, caps), -caps)
    } else {
      w <- pmin(w, caps)
    }
    w
  }

  cap_group_gross <- function(w_row) {
    if (is.null(grp) || is.null(max_per_group)) return(w_row)
    w <- as.numeric(w_row); names(w) <- syms
    split_syms <- split(names(grp), grp)
    for (g in names(split_syms)) {
      idx <- split_syms[[g]]
      gross <- sum(abs(w[idx]), na.rm = TRUE)
      if (is.finite(gross) && gross > max_per_group && gross > 0) {
        w[idx] <- w[idx] * (max_per_group / gross)
      }
    }
    unname(w)
  }

  # row-wise
  M <- as.matrix(DT[, ..syms])
  for (i in seq_len(nrow(M))) {
    w <- M[i, ]
    w <- clamp_symbol(w)
    w <- cap_group_gross(w)

    gross_pos <- sum(pmax(w, 0), na.rm = TRUE)
    if (identical(renormalize, "down") && is.finite(gross_pos) && gross_pos > 1) {
      w <- w / gross_pos
    } else if (identical(renormalize, "both") && is.finite(gross_pos) && gross_pos > 0 && abs(gross_pos - 1) > 1e-12) {
      w <- w / gross_pos
    }

    M[i, ] <- w
  }
  DT[, (syms) := as.data.table(M)]

  if (!is.null(cash_col) && nzchar(cash_col)) {
    if (!(cash_col %in% names(DT))) DT[, (cash_col) := 0]
    gross_pos <- rowSums(pmax(as.matrix(DT[, ..syms]), 0), na.rm = TRUE)
    DT[, (cash_col) := 1 - gross_pos]
  }

  DT[]
}


#' Limit per-date selections to top-K (legacy API)
#'
#' Legacy selector used across examples/vignettes. Works on a WIDE table
#' (`Date` + one column per symbol) and returns a WIDE logical mask with at
#' most `max_positions` TRUE values per row.
#'
#' If `ranking_signal` is supplied, it must have **the same shape** and columns
#' as `selection_df`; the function keeps the top-`K` (largest) scores among
#' the *currently selected* columns in that row. If `ranking_signal` is `NULL`,
#' it falls back to the values in `selection_df` (if numeric), otherwise keeps
#' the first `K` selected symbols in column order (deterministic).
#'
#' @param selection_df data.frame/data.table with columns: `Date`, then one
#'   column per symbol; logical (preferred) or numeric (non-NA/ >0 means selected).
#' @param max_positions positive integer, maximum selections to keep per row.
#' @param ranking_signal optional data.frame/data.table, same dims & columns as
#'   `selection_df`, numeric scores used to rank within the selected set.
#' @param verbose logical; if `TRUE`, prints minor diagnostics. Default `FALSE`.
#'
#' @return A `data.table` with the same columns as `selection_df`, where symbol
#'   columns are logical and at most `max_positions` are TRUE in each row.
#'
#' @examples
#' \donttest{
#'   data(sample_prices_weekly)
#'   mom12 <- PortfolioTesteR::calc_momentum(sample_prices_weekly, 12)  # WIDE numeric
#'   sel10 <- PortfolioTesteR::filter_top_n(mom12, 10)                  # WIDE logical/numeric
#'   # Ensure logical mask
#'   syms <- setdiff(names(sel10), "Date")
#'   sel_mask <- data.table::as.data.table(sel10)
#'   sel_mask[, (syms) := lapply(.SD, function(z) as.logical(as.integer(z))), .SDcols = syms]
#'
#'   # Keep at most 10 per date using momentum as the ranking signal
#'   lim <- limit_positions(selection_df = sel_mask, max_positions = 10L,
#'                          ranking_signal = mom12, verbose = FALSE)
#'   head(lim)
#' }
#'
#' @export
limit_positions <- function(selection_df, max_positions, ranking_signal = NULL, verbose = FALSE) {
  # ---- input checks (self-contained; no reliance on package-internal validators) ----
  if (!is.data.frame(selection_df)) {
    stop("limit_positions : selection_df must be a data.frame or data.table")
  }
  if (!is.numeric(max_positions) || length(max_positions) != 1L || is.na(max_positions) || max_positions <= 0) {
    stop("limit_positions: max_positions must be a positive number")
  }

  DT <- data.table::as.data.table(data.table::copy(selection_df))
  if (!("Date" %in% names(DT))) stop("limit_positions: selection_df must have a 'Date' column")

  syms <- setdiff(names(DT), "Date")
  if (!length(syms)) stop("limit_positions: no symbol columns found in selection_df")

  # Coerce selection to logical mask: numeric>0 -> TRUE; logical kept as-is; others -> FALSE
  DT_mask <- data.table::copy(DT)
  DT_mask[, (syms) := lapply(.SD, function(z) {
    if (is.logical(z)) z
    else if (is.numeric(z)) (is.finite(z) & (z > 0))
    else rep(FALSE, length(z))
  }), .SDcols = syms]

  # Ranking signal (optional)
  use_rank <- !is.null(ranking_signal)
  if (use_rank) {
    RS <- data.table::as.data.table(ranking_signal)
    if (!identical(names(RS), names(DT))) {
      stop("limit_positions: selection_df and ranking_signal must have identical columns (including 'Date')")
    }
    # numeric matrix of scores; non-numeric -> NA
    RS_mat <- as.matrix(data.table::copy(RS[, ..syms]))
    mode(RS_mat) <- "numeric"
  } else {
    # fall back: if selection_df numeric, rank by its numeric values; else NULL
    tmp <- as.matrix(DT[, ..syms])
    if (is.numeric(tmp)) {
      RS_mat <- tmp
      mode(RS_mat) <- "numeric"
      use_rank <- TRUE
    } else {
      RS_mat <- NULL
    }
  }

  M <- as.matrix(DT_mask[, ..syms])  # logical mask
  out <- matrix(FALSE, nrow = nrow(M), ncol = ncol(M))
  colnames(out) <- syms

  K <- as.integer(max_positions)

  for (i in seq_len(nrow(M))) {
    sel_idx <- which(M[i, ])
    if (!length(sel_idx)) {
      # keep none
      next
    }
    if (length(sel_idx) <= K) {
      out[i, sel_idx] <- TRUE
      next
    }
    if (isTRUE(use_rank)) {
      sc <- RS_mat[i, sel_idx]
      sc[!is.finite(sc)] <- -Inf
      top <- sel_idx[order(sc, decreasing = TRUE)][seq_len(K)]
      out[i, top] <- TRUE
    } else {
      # deterministic fallback: first K in column order
      out[i, sel_idx[seq_len(K)]] <- TRUE
    }
  }

  # Build output data.table: same columns, symbol cols as logical
  OUT <- data.table::data.table(Date = DT$Date)
  OUT[, (syms) := as.data.table(as.data.frame(out))]
  data.table::setcolorder(OUT, c("Date", syms))

  if (isTRUE(verbose)) {
    kept_per_row <- rowSums(OUT[, ..syms], na.rm = TRUE)
    cat(sprintf("limit_positions: kept <= %d per row | mean=%.2f max=%d\n",
                K, mean(kept_per_row), max(kept_per_row)))
  }

  OUT[]
}

