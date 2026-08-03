# performance_analytics.R
# Performance analysis with daily monitoring for backtest results
# Calculates enhanced metrics using daily price data

#library(data.table)

#library(zoo)

###############################################################################
# MAIN PERFORMANCE ANALYSIS FUNCTION
###############################################################################


#' Analyze Backtest Performance with Daily Monitoring
#'
#' @description
#' Calculates comprehensive performance metrics using daily price data for
#' enhanced accuracy. Provides risk-adjusted returns, drawdown analysis,
#' and benchmark comparison even when strategy trades at lower frequency.
#'
#' @param backtest_result Result object from run_backtest()
#' @param daily_prices Daily price data including all portfolio symbols
#' @param benchmark_symbol Symbol for benchmark comparison (default: "SPY")
#' @param rf_rate Annual risk-free rate for Sharpe/Sortino (default: 0)
#' @param confidence_level Confidence level for VaR/CVaR (default: 0.95)
#'
#' @return performance_analysis object with metrics and daily tracking
#' @export
#' @examples
#' data("sample_prices_weekly")
#' data("sample_prices_daily")
#'
#' # Use overlapping symbols; cap to 3
#' syms_all <- intersect(names(sample_prices_weekly)[-1], names(sample_prices_daily)[-1])
#' stopifnot(length(syms_all) >= 1)
#' syms <- syms_all[seq_len(min(3L, length(syms_all)))]
#'
#' # Subset weekly (strategy) and daily (monitoring) to the same symbols
#' P <- sample_prices_weekly[, c("Date", syms), with = FALSE]
#' D <- sample_prices_daily[,  c("Date", syms), with = FALSE]
#'
#' # Simple end-to-end example
#' mom <- calc_momentum(P, lookback = 12)
#' sel <- filter_top_n(mom, n = 3)
#' W   <- weight_equally(sel)
#' res <- run_backtest(P, W)
#'
#' # Pick a benchmark that is guaranteed to exist in D
#' perf <- analyze_performance(res, D, benchmark_symbol = syms[1])
#' print(perf)
#' summary(perf)
analyze_performance <- function(backtest_result, daily_prices,
                                benchmark_symbol = "SPY",
                                rf_rate = 0,
                                confidence_level = 0.95) {
  # Analyze backtest performance using daily price monitoring
  #
  # Args:
  #   backtest_result: Result from run_backtest()
  #   daily_prices: Daily price data (must include all symbols + benchmark)
  #   benchmark_symbol: Symbol to use as benchmark (default: SPY)
  #   rf_rate: Annual risk-free rate for Sharpe/Sortino (default: 0)
  #   confidence_level: Confidence level for VaR/CVaR (default: 0.95)
  #
  # Returns:
  #   performance_analysis object with daily tracking and enhanced metrics

  # Input validation
  validate_performance_inputs(backtest_result, daily_prices, benchmark_symbol)

  # Convert to data.table
  setDT(daily_prices)

  # Extract key data from backtest
  positions <- backtest_result$positions
  transactions <- backtest_result$transactions
  strategy_dates <- backtest_result$dates
  initial_capital <- backtest_result$initial_capital

  cat("Analyzing performance with daily monitoring...\n")

  # Step 1: Calculate daily portfolio values
  daily_tracking <- calculate_daily_values(
    positions,
    daily_prices,
    strategy_dates,
    initial_capital,
    backtest_result$cash
  )

  # Step 2: Calculate daily returns
  daily_returns <- c(0, diff(daily_tracking$portfolio_values) /
                       head(daily_tracking$portfolio_values, -1))

  # Step 3: Calculate enhanced metrics from daily data
  metrics <- calculate_enhanced_metrics(
    daily_values = daily_tracking$portfolio_values,
    daily_returns = daily_returns,
    rf_rate = rf_rate,
    confidence_level = confidence_level
  )

  # Step 4: Benchmark comparison
  benchmark_analysis <- NULL
  if (benchmark_symbol %in% names(daily_prices)) {
    # FIXED: Pass full benchmark data with dates instead of just price vector
    benchmark_data <- daily_prices[, c("Date", benchmark_symbol), with = FALSE]

    benchmark_analysis <- analyze_vs_benchmark(
      portfolio_returns = daily_returns,
      benchmark_prices = benchmark_data,  # Changed: pass full dataframe
      dates = daily_tracking$dates,
      benchmark_symbol = benchmark_symbol
    )
  } else {
    warning(sprintf("Benchmark '%s' not found in daily_prices. Skipping benchmark analysis.",
                    benchmark_symbol))
  }

  # Step 5: Period-based analysis
  period_analysis <- analyze_by_period(
    dates = daily_tracking$dates,
    returns = daily_returns,
    values = daily_tracking$portfolio_values
  )

  # Step 6: Risk analysis
  risk_analysis <- calculate_risk_metrics(
    returns = daily_returns,
    values = daily_tracking$portfolio_values,
    confidence_level = confidence_level
  )

  # Create performance analysis object
  result <- list(
    # Original backtest info
    strategy_name = backtest_result$name,
    strategy_frequency = get_data_frequency(strategy_dates),
    monitoring_frequency = "daily",

    # Daily tracking
    daily = list(
      dates = daily_tracking$dates,
      values = daily_tracking$portfolio_values,
      returns = daily_returns,
      positions_value = daily_tracking$positions_value,
      cash = daily_tracking$cash
    ),

    # Enhanced metrics
    metrics = metrics,

    # Risk analysis
    risk = risk_analysis,

    # Benchmark comparison
    benchmark = benchmark_analysis,

    # Period analysis
    periods = period_analysis,

    # Original data (for reference)
    original_result = backtest_result,
    daily_prices = daily_prices
  )

  class(result) <- c("performance_analysis", "list")

  return(result)
}

###############################################################################
# DAILY VALUE CALCULATION
###############################################################################
#' Daily equity curve from positions and daily prices
#'
#' Carries portfolio positions (from a weekly or lower-frequency backtest)
#' forward to daily dates, multiplies by daily prices, and combines with cash
#' to produce a daily portfolio value series for monitoring and analytics.
#'
#' @param positions A `data.frame`/`data.table` of portfolio positions with
#'   columns `Date` + symbols. Values should be the backtest's **position
#'   inventory** per symbol at each rebalance date (typically shares or notional
#'   units consistent with your backtest's accounting).
#' @param daily_prices A `data.frame`/`data.table` of **daily** prices with
#'   columns `Date` + the same symbol set present in `positions` (at least the
#'   intersection).
#' @param strategy_dates A `Date` vector of the backtest's decision/rebalance
#'   calendar (e.g., `backtest_result$dates`).
#' @param initial_capital Numeric scalar. Starting cash used for days **before**
#'   the first position exists (typically `backtest_result$initial_capital`).
#' @param cash_values Optional numeric vector of cash balances at the strategy
#'   dates (e.g., `backtest_result$cash`). If `NULL`, leading days are treated
#'   as all-cash (= `initial_capital`) and post-rebalance cash defaults to 0.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{dates} Daily dates within the strategy span.
#'   \item \code{portfolio_values} Daily total portfolio value (positions + cash).
#'   \item \code{positions_value} Daily mark-to-market of positions only.
#'   \item \code{cash} Daily carried cash series.
#' }
#'
#' @examples
#' \donttest{
#'   # Minimal end-to-end example using bundled data and a simple weekly backtest
#'   library(PortfolioTesteR)
#'   data(sample_prices_weekly); data(sample_prices_daily)
#'
#'   # Build a tiny strategy: momentum -> top-3 -> equal weights
#'   mom <- calc_momentum(sample_prices_weekly, lookback = 12)
#'   sel <- filter_top_n(mom, n = 3)
#'   W   <- weight_equally(sel)
#'   bt  <- run_backtest(sample_prices_weekly, W, name = "Demo")
#'
#'   # Compute daily monitoring values from positions + cash
#'   vals <- calculate_daily_values(
#'     positions       = bt$positions,
#'     daily_prices    = sample_prices_daily,
#'     strategy_dates  = bt$dates,
#'     initial_capital = bt$initial_capital,
#'     cash_values     = bt$cash
#'   )
#'
#'   # Quick sanity checks
#'   head(vals$dates)
#'   head(vals$portfolio_values)
#' }
#'
#' @export
calculate_daily_values <- function(positions,
                                   daily_prices,
                                   strategy_dates,
                                   initial_capital,
                                   cash_values) {

  positions    <- data.table::as.data.table(positions)
  daily_prices <- data.table::as.data.table(daily_prices)
  data.table::setorder(positions, Date)
  data.table::setorder(daily_prices, Date)

  # Symbols taken from positions (backtest output)
  syms <- setdiff(names(positions), "Date")

  # Restrict daily window to strategy span and keep only needed symbols
  dmin <- min(strategy_dates)
  dmax <- max(strategy_dates)
  dp   <- daily_prices[Date >= dmin & Date <= dmax, c("Date", syms), with = FALSE]

  # Carry forward positions to each daily date (roll join on Date)
  pos_cf <- positions[list(Date = dp$Date), on = "Date", roll = TRUE]

  # If there are no positions at all, treat as all-cash every day
  if (nrow(pos_cf) == 0L) {
    pos_val     <- rep(0, nrow(dp))
    cash_series <- rep(initial_capital, nrow(dp))
    total       <- cash_series
    return(list(
      dates            = dp$Date,
      portfolio_values = total,
      positions_value  = pos_val,
      cash             = cash_series
    ))
  }

  # Ensure all symbol columns exist; fill missing with 0 (no position)
  missing_syms <- setdiff(syms, names(pos_cf))
  for (s in missing_syms) pos_cf[[s]] <- NA_real_
  pos_cf[, (syms) := lapply(.SD, function(x) data.table::fifelse(is.na(x), 0, x)),
         .SDcols = syms]

  # Mark-to-market of positions
  P <- as.matrix(dp[, ..syms])      # daily prices
  W <- as.matrix(pos_cf[, ..syms])  # carried-forward positions (shares/notional)
  pos_val <- rowSums(W * P, na.rm = TRUE)

  # Carry cash to daily dates
  if (is.null(cash_values)) {
    cash_series <- rep(0, nrow(dp))
    first_pos_date <- positions$Date[1L]
    pre_mask <- dp$Date < first_pos_date
    if (any(pre_mask)) cash_series[pre_mask] <- initial_capital
  } else {
    cash_dt <- data.table::data.table(Date = as.Date(strategy_dates),
                                      cash = as.numeric(cash_values))
    data.table::setkey(cash_dt, Date)
    cash_cf <- cash_dt[list(Date = dp$Date), roll = TRUE]
    cash_series <- cash_cf$cash
    if (length(cash_series) && is.na(cash_series[1L])) cash_series[1L] <- initial_capital
    cash_series[is.na(cash_series)] <- 0
  }

  total <- pos_val + cash_series

  list(
    dates            = dp$Date,
    portfolio_values = total,
    positions_value  = pos_val,
    cash             = cash_series
  )
}



###############################################################################
# ENHANCED METRICS CALCULATION
###############################################################################
#' Calculate Enhanced Performance Metrics
#'
#' @description
#' Computes comprehensive risk and return metrics from daily data including
#' Sharpe, Sortino, Calmar ratios, VaR, CVaR, and tail risk measures.
#'
#' @param daily_values Daily portfolio values
#' @param daily_returns Daily return series
#' @param rf_rate Risk-free rate
#' @param confidence_level VaR/CVaR confidence level
#'
#' @return List of performance metrics
#' @keywords internal
calculate_enhanced_metrics <- function(daily_values, daily_returns,
                                       rf_rate = 0, confidence_level = 0.95) {
  # Calculate comprehensive metrics from daily data

  # Annualization factor (252 trading days)
  annual_factor <- 252

  # Basic return metrics
  total_return <- (tail(daily_values, 1) / head(daily_values, 1)) - 1
  n_days <- length(daily_values)
  n_years <- n_days / annual_factor
  annualized_return <- (1 + total_return)^(1/n_years) - 1

  # Volatility metrics
  daily_vol <- sd(daily_returns, na.rm = TRUE)
  annual_vol <- daily_vol * sqrt(annual_factor)

  # Downside volatility (for Sortino)
  downside_returns <- daily_returns[daily_returns < rf_rate/annual_factor]
  downside_vol <- sd(downside_returns, na.rm = TRUE) * sqrt(annual_factor)

  # Risk-adjusted returns
  sharpe_ratio <- (annualized_return - rf_rate) / annual_vol
  sortino_ratio <- (annualized_return - rf_rate) / downside_vol

  # Drawdown analysis
  cumulative <- cumprod(1 + daily_returns)
  running_max <- cummax(cumulative)
  drawdown <- (cumulative - running_max) / running_max
  max_drawdown <- min(drawdown)

  # Calmar ratio (return / max drawdown)
  calmar_ratio <- annualized_return / abs(max_drawdown)

  # Value at Risk (VaR) and Conditional VaR (CVaR)
  var_daily <- quantile(daily_returns, 1 - confidence_level)
  cvar_daily <- mean(daily_returns[daily_returns <= var_daily])
  var_annual <- var_daily * sqrt(annual_factor)
  cvar_annual <- cvar_daily * sqrt(annual_factor)

  # Win rate and average win/loss
  winning_days <- sum(daily_returns > 0)
  total_days <- length(daily_returns)
  daily_win_rate <- winning_days / total_days

  avg_win <- mean(daily_returns[daily_returns > 0])
  avg_loss <- mean(daily_returns[daily_returns < 0])
  profit_factor <- abs(sum(daily_returns[daily_returns > 0]) /
                         sum(daily_returns[daily_returns < 0]))

  # Compile metrics
  metrics <- list(
    # Returns
    total_return = total_return,
    annualized_return = annualized_return,

    # Volatility
    daily_volatility = daily_vol,
    annual_volatility = annual_vol,
    downside_volatility = downside_vol,

    # Risk-adjusted
    sharpe_ratio = sharpe_ratio,
    sortino_ratio = sortino_ratio,
    calmar_ratio = calmar_ratio,

    # Drawdown
    max_drawdown = max_drawdown,
    avg_drawdown = mean(drawdown[drawdown < 0]),

    # Risk metrics
    var_95_daily = var_daily,
    cvar_95_daily = cvar_daily,
    var_95_annual = var_annual,
    cvar_95_annual = cvar_annual,

    # Win/loss
    daily_win_rate = daily_win_rate,
    avg_daily_win = avg_win,
    avg_daily_loss = avg_loss,
    profit_factor = profit_factor
  )

  return(metrics)
}

###############################################################################
# BENCHMARK COMPARISON
###############################################################################

#' Benchmark-relative performance statistics
#'
#' Computes standard benchmark-relative metrics (e.g., correlation, beta, alpha, tracking error,
#' information ratio) by aligning portfolio returns with benchmark returns derived from prices.
#'
#' @param portfolio_returns A numeric vector of portfolio simple returns aligned to \code{dates}.
#' @param benchmark_prices A data frame (Date + symbols) of adjusted benchmark prices at the
#'   same cadence as \code{dates}.
#' @param dates A vector of \code{Date} values used to align \code{portfolio_returns}
#'   with \code{benchmark_prices}.
#' @param benchmark_symbol Character scalar giving the column name (symbol) in \code{benchmark_prices}
#'   to use as the benchmark.
#'
#' @return A list or data frame with benchmark-relative statistics according to the package's
#'   conventions, including correlation, beta, alpha, tracking error, and information ratio.
#'
#' @examples
#' \donttest{
#'   data(sample_prices_weekly)
#'   mom12 <- PortfolioTesteR::calc_momentum(sample_prices_weekly, lookback = 12)
#'   sel10 <- PortfolioTesteR::filter_top_n(mom12, n = 5)
#'   w_eq  <- PortfolioTesteR::weight_equally(sel10)
#'   pr    <- PortfolioTesteR::portfolio_returns(w_eq, sample_prices_weekly)
#'
#'   # Use SPY as the benchmark
#'   bench <- sample_prices_weekly[, c("Date", "SPY")]
#'   res <- analyze_vs_benchmark(
#'     pr$portfolio_return,
#'     bench,
#'     dates = pr$Date,
#'     benchmark_symbol = "SPY"
#'   )
#'   res
#' }
#'
#' @export
analyze_vs_benchmark <- function(portfolio_returns, benchmark_prices, dates, benchmark_symbol = "SPY") {
  # Compare portfolio performance against benchmark

  # Extract price column from dataframe
  if (is.data.frame(benchmark_prices) || is.data.table(benchmark_prices)) {
    benchmark_dates <- benchmark_prices$Date
    benchmark_price_vector <- benchmark_prices[[benchmark_symbol]]
  } else {
    # Fallback for backward compatibility
    benchmark_price_vector <- benchmark_prices
    benchmark_dates <- dates  # Assume aligned
  }

  # CRITICAL FIX: Align benchmark data to portfolio dates
  # Create a data.table for easy date joining
  portfolio_dt <- data.table(
    date = dates,
    portfolio_return = portfolio_returns
  )

  benchmark_dt <- data.table(
    date = benchmark_dates,
    benchmark_price = benchmark_price_vector
  )

  # Calculate benchmark returns
  benchmark_dt[, benchmark_return := c(0, diff(benchmark_price) / head(benchmark_price, -1))]

  # Join on dates to ensure alignment
  aligned_data <- portfolio_dt[benchmark_dt, on = "date", nomatch = 0]

  # Extract aligned returns
  portfolio_returns_clean <- aligned_data$portfolio_return
  benchmark_returns_clean <- aligned_data$benchmark_return

  # Remove any NA values
  valid_idx <- !is.na(portfolio_returns_clean) & !is.na(benchmark_returns_clean)
  portfolio_returns_clean <- portfolio_returns_clean[valid_idx]
  benchmark_returns_clean <- benchmark_returns_clean[valid_idx]

  # Check if we have enough data
  if(length(portfolio_returns_clean) < 10) {
    warning("Not enough valid observations for reliable benchmark comparison")
    return(NULL)
  }

  # Now calculate correlation and beta with properly aligned data
  correlation <- cor(portfolio_returns_clean, benchmark_returns_clean)

  # Beta (portfolio sensitivity to benchmark)
  covariance <- cov(portfolio_returns_clean, benchmark_returns_clean)
  benchmark_variance <- var(benchmark_returns_clean)
  beta <- covariance / benchmark_variance

  # Alpha (excess return)
  # Using CAPM: Portfolio Return = Alpha + Beta * Benchmark Return
  portfolio_annual <- mean(portfolio_returns_clean) * 252
  benchmark_annual <- mean(benchmark_returns_clean) * 252
  alpha <- portfolio_annual - beta * benchmark_annual

  # Tracking error (standard deviation of excess returns)
  excess_returns <- portfolio_returns_clean - benchmark_returns_clean
  tracking_error <- sd(excess_returns) * sqrt(252)

  # Information ratio (active return / tracking error)
  active_return <- portfolio_annual - benchmark_annual
  information_ratio <- active_return / tracking_error

  # Relative performance
  portfolio_cumulative <- prod(1 + portfolio_returns_clean) - 1
  benchmark_cumulative <- prod(1 + benchmark_returns_clean) - 1
  relative_performance <- portfolio_cumulative - benchmark_cumulative

  return(list(
    symbol = benchmark_symbol,
    correlation = correlation,
    beta = beta,
    alpha = alpha,
    tracking_error = tracking_error,
    information_ratio = information_ratio,
    portfolio_total_return = portfolio_cumulative,
    benchmark_total_return = benchmark_cumulative,
    relative_performance = relative_performance,
    excess_returns = excess_returns,
    n_observations = length(portfolio_returns_clean)
  ))
}

#' Period-level summary statistics
#'
#' Aggregates portfolio results by calendar period and computes standard statistics
#' for each period. Provide at least one of `returns` or `values`.
#'
#' @param dates Date vector aligned to `returns` / `values`.
#' @param returns Numeric simple returns aligned to `dates` (optional).
#' @param values Numeric equity values aligned to `dates` (optional).
#' @param period "monthly", "quarterly", or "yearly".
#' @param na_rm Logical; remove NAs inside per-period aggregations.
#' @return data.frame with period keys and columns: ret, start_value, end_value, n_obs.
#' @examples
#' \donttest{
#'   data(sample_prices_weekly)
#'   mom12 <- PortfolioTesteR::calc_momentum(sample_prices_weekly, lookback = 12)
#'   sel5  <- PortfolioTesteR::filter_top_n(mom12, n = 5)
#'   w_eq  <- PortfolioTesteR::weight_equally(sel5)
#'   pr    <- PortfolioTesteR::portfolio_returns(w_eq, sample_prices_weekly)
#'   val   <- 1e5 * cumprod(1 + pr$portfolio_return)
#'   out   <- analyze_by_period(
#'     dates   = pr$Date,
#'     returns = pr$portfolio_return,
#'     values  = val,
#'     period  = "monthly"
#'   )
#'   head(out)
#' }
#' @export
analyze_by_period <- function(dates,
                              returns = NULL,
                              values  = NULL,
                              period  = c("monthly","quarterly","yearly"),
                              na_rm   = TRUE) {
  period <- match.arg(period)
  if (length(dates) == 0L) stop("analyze_by_period: `dates` is empty.")
  if (is.null(returns) && is.null(values)) {
    stop("analyze_by_period: provide at least one of `returns` or `values`.")
  }

  dt <- data.table::data.table(
    Date    = as.Date(dates),
    ret_col = if (!is.null(returns)) as.numeric(returns) else NA_real_,
    val_col = if (!is.null(values))  as.numeric(values)  else NA_real_
  )
  has_ret <- !is.null(returns)
  has_val <- !is.null(values)

  comp_ret <- function(x) {
    x <- as.numeric(x); if (na_rm) x <- x[!is.na(x)]
    if (!length(x)) NA_real_ else prod(1 + x) - 1
  }
  first_num <- function(x) {
    x <- as.numeric(x); i <- which(!is.na(x))[1L]; if (is.na(i)) NA_real_ else x[i]
  }
  last_num <- function(x) {
    x <- as.numeric(x); idx <- which(!is.na(x)); if (!length(idx)) NA_real_ else x[idx[length(idx)]]
  }

  if (identical(period, "monthly")) {
    res <- dt[
      , {
        r <- .SD[["ret_col"]]; v <- .SD[["val_col"]]
        list(
          ret         = if (has_ret) comp_ret(r) else NA_real_,
          start_value = if (has_val) first_num(v) else NA_real_,
          end_value   = if (has_val) last_num(v)  else NA_real_,
          n_obs       = .N
        )
      },
      by = list(Year  = as.integer(format(Date, "%Y")),
                Month = as.integer(format(Date, "%m"))),
      .SDcols = c("ret_col","val_col")
    ]
  } else if (identical(period, "quarterly")) {
    res <- dt[
      , {
        r <- .SD[["ret_col"]]; v <- .SD[["val_col"]]
        list(
          ret         = if (has_ret) comp_ret(r) else NA_real_,
          start_value = if (has_val) first_num(v) else NA_real_,
          end_value   = if (has_val) last_num(v)  else NA_real_,
          n_obs       = .N
        )
      },
      by = list(Year    = as.integer(format(Date, "%Y")),
                Quarter = ((as.integer(format(Date, "%m")) - 1L) %/% 3L) + 1L),
      .SDcols = c("ret_col","val_col")
    ]
  } else { # yearly
    res <- dt[
      , {
        r <- .SD[["ret_col"]]; v <- .SD[["val_col"]]
        list(
          ret         = if (has_ret) comp_ret(r) else NA_real_,
          start_value = if (has_val) first_num(v) else NA_real_,
          end_value   = if (has_val) last_num(v)  else NA_real_,
          n_obs       = .N
        )
      },
      by = list(Year = as.integer(format(Date, "%Y"))),
      .SDcols = c("ret_col","val_col")
    ]
  }

  out <- as.data.frame(res)
  rownames(out) <- NULL
  out
}


#######################
###############################################################################
# RISK METRICS
###############################################################################

calculate_risk_metrics <- function(returns, values, confidence_level = 0.95) {
  # Calculate comprehensive risk metrics

  # Drawdown analysis
  drawdowns <- calculate_drawdown_series(values)
  dd_info <- analyze_drawdowns(drawdowns, returns)

  # Volatility analysis
  rolling_vol_30d <- zoo::rollapply(returns, 30, sd, fill = NA, align = "right")
  rolling_vol_90d <- zoo::rollapply(returns, 90, sd, fill = NA, align = "right")

  # Tail risk
  left_tail <- quantile(returns, probs = c(0.01, 0.05, 0.10), na.rm = TRUE)
  right_tail <- quantile(returns, probs = c(0.90, 0.95, 0.99), na.rm = TRUE)

  # Skewness and kurtosis
  skewness <- sum((returns - mean(returns))^3) / (length(returns) * sd(returns)^3)
  kurtosis <- sum((returns - mean(returns))^4) / (length(returns) * sd(returns)^4) - 3

  return(list(
    drawdown_info = dd_info,
    volatility = list(
      current = tail(rolling_vol_30d[!is.na(rolling_vol_30d)], 1),
      vol_30d_avg = mean(rolling_vol_30d, na.rm = TRUE),
      vol_90d_avg = mean(rolling_vol_90d, na.rm = TRUE),
      vol_percentile_95 = quantile(rolling_vol_30d, 0.95, na.rm = TRUE)
    ),
    tail_risk = list(
      left_tail = left_tail,
      right_tail = right_tail,
      skewness = skewness,
      kurtosis = kurtosis
    )
  ))
}

###############################################################################
# UTILITY FUNCTIONS
###############################################################################
#' Calculate Drawdown Time Series
#'
#' @description
#' Computes drawdown series from portfolio values.
#'
#' @param values Numeric vector of portfolio values
#'
#' @return Numeric vector of drawdowns (as negative percentages)
#' @export
#' @examples
#' data("sample_prices_weekly")
#' momentum <- calc_momentum(sample_prices_weekly, lookback = 12)
#' sel <- filter_top_n(momentum, n = 10)
#' W   <- weight_equally(sel)
#' res <- run_backtest(sample_prices_weekly, W)
#' dd_series <- calculate_drawdown_series(res$portfolio_values)
#' dd_stats  <- analyze_drawdowns(dd_series, res$returns)
calculate_drawdown_series <- function(values) {
  cummax_values <- cummax(values)
  drawdowns <- (values - cummax_values) / cummax_values
  return(drawdowns)
}

#' Analyze Drawdown Characteristics
#'
#' @description
#' Detailed analysis of drawdown periods including depth, duration, and recovery.
#'
#' @param drawdowns Drawdown series (negative values)
#' @param returns Return series for additional metrics
#'
#' @return List with drawdown statistics
#' @export
#' @examples
#' data("sample_prices_weekly")
#' momentum <- calc_momentum(sample_prices_weekly, lookback = 12)
#' selected <- filter_top_n(momentum, n = 10)
#' weights <- weight_equally(selected)
#' result <- run_backtest(sample_prices_weekly, weights)
#' dd_analysis <- analyze_drawdowns(result$portfolio_value, result$dates)
analyze_drawdowns <- function(drawdowns, returns) {
  # Find drawdown periods
  in_drawdown <- drawdowns < 0

  # Identify start and end of each drawdown period
  starts <- which(diff(c(FALSE, in_drawdown)) == 1)
  ends <- which(diff(c(in_drawdown, FALSE)) == -1)

  if (length(starts) == 0) {
    return(list(
      max_drawdown = 0,
      avg_drawdown = 0,
      max_duration = 0,
      avg_duration = 0,
      current_drawdown = 0,
      n_drawdowns = 0
    ))
  }

  # Calculate drawdown statistics
  drawdown_depths <- numeric(length(starts))
  drawdown_durations <- ends - starts + 1

  for (i in seq_along(starts)) {
    drawdown_depths[i] <- min(drawdowns[starts[i]:ends[i]])
  }

  return(list(
    max_drawdown = min(drawdowns),
    avg_drawdown = mean(drawdown_depths),
    max_duration = max(drawdown_durations),
    avg_duration = mean(drawdown_durations),
    current_drawdown = tail(drawdowns, 1),
    n_drawdowns = length(starts),
    recovery_time = calculate_recovery_time(drawdowns)
  ))
}

calculate_recovery_time <- function(drawdowns) {
  # Time to recover from max drawdown
  max_dd_idx <- which.min(drawdowns)
  if (max_dd_idx == length(drawdowns)) {
    return(NA)  # Still in max drawdown
  }

  # Find recovery point
  recovery_idx <- which(drawdowns[(max_dd_idx+1):length(drawdowns)] >= 0)[1]
  if (is.na(recovery_idx)) {
    return(NA)  # Not recovered yet
  }

  return(recovery_idx)
}



#' Detect Data Frequency from Dates
#'
#' @description
#' Automatically detects whether data is daily, weekly, monthly, or
#' quarterly based on date spacing.
#'
#' @param dates Vector of Date objects
#'
#' @return Character string: "daily", "weekly", "monthly", or "quarterly"
#' @export
#' @examples
#' data("sample_prices_weekly")
#' freq <- get_data_frequency(sample_prices_weekly$Date)
get_data_frequency <- function(dates) {
  # Detect data frequency from dates
  date_diffs <- as.numeric(diff(dates))
  median_diff <- median(date_diffs)

  if (median_diff <= 1) return("daily")
  if (median_diff <= 7) return("weekly")
  if (median_diff <= 31) return("monthly")
  return("quarterly")
}

###############################################################################
# VALIDATION
###############################################################################

#' Validate Performance Analysis Inputs
#'
#' @description
#' data("sample_prices_weekly")
#' Checks that backtest result and daily prices are properly formatted
#' with matching symbols and appropriate date coverage.
#'
#' @param backtest_result Backtest result object
#' data("sample_prices_weekly")
#' @param daily_prices Daily price data
#' @param benchmark_symbol Benchmark symbol
#'
#' @return TRUE if valid, stops with error if not
#' @keywords internal
validate_performance_inputs <- function(backtest_result, daily_prices, benchmark_symbol) {
  # Validate inputs

  if (!inherits(backtest_result, "backtest_result")) {
    stop("backtest_result must be output from run_backtest()")
  }

  if (!is.data.frame(daily_prices)) {
    stop("daily_prices must be a data.frame or data.table")
  }

  # Check if daily frequency
  dates <- daily_prices$Date
  date_diffs <- as.numeric(diff(dates))
  median_diff <- median(date_diffs)

  if (median_diff > 5) {
    stop(sprintf("daily_prices appears to be %s frequency, not daily.
                 Daily data is required for performance analysis.",
                 ifelse(median_diff > 20, "monthly", "weekly")))
  }

  # Check if all portfolio symbols exist in daily data
  position_symbols <- setdiff(names(backtest_result$positions), "Date")
  missing_symbols <- position_symbols[!position_symbols %in% names(daily_prices)]

  if (length(missing_symbols) > 0) {
    stop(sprintf("The following symbols from portfolio are missing in daily_prices: %s",
                 paste(missing_symbols, collapse = ", ")))
  }

  # Check date coverage
  strategy_start <- min(backtest_result$dates)
  strategy_end <- max(backtest_result$dates)
  daily_start <- min(daily_prices$Date)
  daily_end <- max(daily_prices$Date)

  if (daily_start > strategy_start || daily_end < strategy_end) {
    stop(sprintf("daily_prices date range (%s to %s) does not fully cover strategy period (%s to %s)",
                 daily_start, daily_end, strategy_start, strategy_end))
  }

  return(TRUE)
}

###############################################################################
# S3 METHODS
###############################################################################
#' Print Performance Analysis Results
#'
#' @description
#' S3 method for printing performance analysis with key metrics including
#' risk-adjusted returns, drawdown statistics, and benchmark comparison.
#'
#' @param x performance_analysis object
#' @param ... Additional arguments (unused)
#'
#' @return Invisible copy of x
#' @export
#' @examples
#' data("sample_prices_weekly")
#' data("sample_prices_daily")
#' syms_all <- intersect(names(sample_prices_weekly)[-1], names(sample_prices_daily)[-1])
#' syms <- syms_all[seq_len(min(3L, length(syms_all)))]
#' P <- sample_prices_weekly[, c("Date", syms), with = FALSE]
#' D <- sample_prices_daily[,  c("Date", syms), with = FALSE]
#' mom <- calc_momentum(P, lookback = 12)
#' sel <- filter_top_n(mom, n = 3)
#' W   <- weight_equally(sel)
#' res <- run_backtest(P, W)
#' perf <- analyze_performance(res, D, benchmark_symbol = syms[1])
#' print(perf)  # or just: perf
print.performance_analysis <- function(x, ...) {
  cat("Performance Analysis:", x$strategy_name, "\n")
  cat("=========================================\n")
  cat(sprintf("Strategy Frequency: %s\n", x$strategy_frequency))
  cat(sprintf("Monitoring Frequency: %s (enhanced metrics)\n", x$monitoring_frequency))
  cat(sprintf("Analysis Period: %s to %s\n",
              min(x$daily$dates), max(x$daily$dates)))

  cat("\nRisk-Adjusted Performance:\n")
  cat(sprintf("  Total Return: %.2f%%\n", x$metrics$total_return * 100))
  cat(sprintf("  Annualized Return: %.2f%%\n", x$metrics$annualized_return * 100))
  cat(sprintf("  Sharpe Ratio: %.2f\n", x$metrics$sharpe_ratio))
  cat(sprintf("  Sortino Ratio: %.2f\n", x$metrics$sortino_ratio))
  cat(sprintf("  Calmar Ratio: %.2f\n", x$metrics$calmar_ratio))

  cat("\nRisk Metrics:\n")
  cat(sprintf("  Max Drawdown: %.2f%%\n", x$metrics$max_drawdown * 100))
  cat(sprintf("  Annual Volatility: %.2f%%\n", x$metrics$annual_volatility * 100))
  cat(sprintf("  Daily 95%% VaR: %.2f%%\n", x$metrics$var_95_daily * 100))
  cat(sprintf("  Daily Win Rate: %.1f%%\n", x$metrics$daily_win_rate * 100))

  if (!is.null(x$benchmark)) {
    cat("\nBenchmark Comparison:\n")
    cat(sprintf("  Beta: %.2f\n", x$benchmark$beta))
    cat(sprintf("  Alpha: %.2f%%\n", x$benchmark$alpha * 100))
    cat(sprintf("  Correlation: %.2f\n", x$benchmark$correlation))
    cat(sprintf("  Tracking Error: %.2f%%\n", x$benchmark$tracking_error * 100))
  }

  invisible(x)
}

#' Plot Performance Analysis Results
#'
#' @description
#' S3 method for visualizing performance metrics. Supports multiple plot
#' types including summary dashboard, return distributions, risk evolution,
#' and rolling statistics.
#'
#' @param x performance_analysis object
#' @param type Plot type: "summary", "returns", "risk", "drawdown"
#' @param ... Additional plotting parameters
#'
#' @return NULL (creates plot)
#' @export
#' @examples
#' data("sample_prices_weekly")
#' data("sample_prices_daily")
#' syms_all <- intersect(names(sample_prices_weekly)[-1], names(sample_prices_daily)[-1])
#' syms <- syms_all[seq_len(min(3L, length(syms_all)))]
#' P <- sample_prices_weekly[, c("Date", syms), with = FALSE]
#' D <- sample_prices_daily[,  c("Date", syms), with = FALSE]
#' mom <- calc_momentum(P, lookback = 12)
#' sel <- filter_top_n(mom, n = 3)
#' W   <- weight_equally(sel)
#' res <- run_backtest(P, W)
#' perf <- analyze_performance(res, D, benchmark_symbol = syms[1])
#' if (interactive()) {
#'   plot(perf, type = "summary")
#' }
plot.performance_analysis <- function(x, type = "summary", ...) {
  # Save and restore only mfrow
  old_mfrow <- par("mfrow")
  on.exit(par(mfrow = old_mfrow))

  # Create various plots for performance analysis
  if (type == "summary") {
    # 4-panel summary plot
    par(mfrow = c(2, 2))

    # Panel 1: Portfolio value
    plot(x$daily$dates, x$daily$values, type = "l",
         col = "darkblue", lwd = 2,
         main = "Daily Portfolio Value",
         xlab = "Date", ylab = "Value ($)")
    grid()

    # Panel 2: Drawdowns
    dd <- calculate_drawdown_series(x$daily$values)
    plot(x$daily$dates, dd * 100, type = "l",
         col = "darkred", lwd = 2,
         main = "Drawdowns",
         xlab = "Date", ylab = "Drawdown (%)")
    abline(h = 0, lty = 2)
    grid()

    # Panel 3: Rolling volatility
    if (length(x$daily$returns) >= 30) {
      vol_30d <- numeric(length(x$daily$returns))
      for (i in 30:length(x$daily$returns)) {
        vol_30d[i] <- sd(x$daily$returns[(i-29):i]) * sqrt(252) * 100
      }
      vol_30d[1:29] <- NA

      plot(x$daily$dates, vol_30d, type = "l",
           col = "darkgreen", lwd = 2,
           main = "30-Day Rolling Volatility",
           xlab = "Date", ylab = "Annualized Vol (%)")
      grid()
    }

    # Panel 4: Cumulative returns
    cum_ret <- cumprod(1 + x$daily$returns) - 1
    plot(x$daily$dates, cum_ret * 100, type = "l",
         col = "darkblue", lwd = 2,
         main = "Cumulative Returns",
         xlab = "Date", ylab = "Return (%)")

    if (!is.null(x$benchmark)) {
      # Add benchmark if available
      bench_cum <- cumprod(1 + x$benchmark$excess_returns + x$daily$returns) - 1
      lines(x$daily$dates[1:length(bench_cum)], bench_cum * 100,
            col = "gray50", lwd = 2, lty = 2)
      legend("topleft", c("Portfolio", "Benchmark"),
             col = c("darkblue", "gray50"), lwd = 2, lty = c(1, 2))
    }
    grid()

  } else if (type == "returns") {
    # Return distribution plots
    par(mfrow = c(2, 2))

    # Histogram
    hist(x$daily$returns * 100, breaks = 50,
         main = "Daily Return Distribution",
         xlab = "Daily Return (%)", col = "lightblue")
    abline(v = mean(x$daily$returns) * 100, col = "red", lwd = 2)

    # Q-Q plot
    qqnorm(x$daily$returns)
    qqline(x$daily$returns, col = "red")

    # ACF
    acf(x$daily$returns, main = "Return Autocorrelation")

    # Monthly returns - FIXED column name
    if (!is.null(x$periods$monthly)) {
      barplot(x$periods$monthly$ret * 100,  # Fixed: was 'return', now 'ret'
              names.arg = x$periods$monthly$year_month,
              las = 2, cex.names = 0.6,
              main = "Monthly Returns",
              ylab = "Return (%)",
              col = ifelse(x$periods$monthly$ret > 0, "darkgreen", "darkred"))
      abline(h = 0)
    }

  } else if (type == "risk") {
    # Risk-focused plots
    par(mfrow = c(2, 2))

    # VaR visualization
    returns_sorted <- sort(x$daily$returns)
    n <- length(returns_sorted)
    plot(1:n / n, returns_sorted * 100, type = "l",
         main = "Return Distribution (CDF)",
         xlab = "Percentile", ylab = "Daily Return (%)")
    abline(h = x$metrics$var_95_daily * 100, col = "red", lty = 2)
    abline(v = 0.05, col = "red", lty = 2)
    text(0.05, x$metrics$var_95_daily * 100 + 0.5, "95% VaR", col = "red")
    grid()

    # Drawdown periods
    dd <- x$risk$drawdown_info
    plot(x$daily$dates, calculate_drawdown_series(x$daily$values) * 100,
         type = "l", col = "darkred",
         main = sprintf("Drawdown Analysis (Max: %.1f%%)", dd$max_drawdown * 100),
         xlab = "Date", ylab = "Drawdown (%)")
    abline(h = dd$avg_drawdown * 100, lty = 2, col = "red")
    text(x$daily$dates[1], dd$avg_drawdown * 100 - 1, "Avg DD", col = "red")
    grid()

    # Rolling Sharpe
    if (length(x$daily$returns) >= 252) {
      rolling_sharpe <- numeric(length(x$daily$returns))
      for (i in 252:length(x$daily$returns)) {
        period_returns <- x$daily$returns[(i-251):i]
        rolling_sharpe[i] <- mean(period_returns) / sd(period_returns) * sqrt(252)
      }
      rolling_sharpe[1:251] <- NA

      plot(x$daily$dates, rolling_sharpe, type = "l",
           col = "purple", lwd = 2,
           main = "Rolling 1-Year Sharpe Ratio",
           xlab = "Date", ylab = "Sharpe Ratio")
      abline(h = x$metrics$sharpe_ratio, lty = 2)
      grid()
    }

    # Beta over time (if benchmark available)
    if (!is.null(x$benchmark) && length(x$daily$returns) >= 60) {
      rolling_beta <- numeric(length(x$daily$returns))
      bench_returns <- c(0, diff(x$daily_prices$SPY) / head(x$daily_prices$SPY, -1))

      for (i in 60:min(length(x$daily$returns), length(bench_returns))) {
        period_port <- x$daily$returns[(i-59):i]
        period_bench <- bench_returns[(i-59):i]
        if (length(period_port) == length(period_bench)) {
          cov_val <- cov(period_port, period_bench, use = "complete.obs")
          var_bench <- var(period_bench, na.rm = TRUE)
          if (var_bench > 0) {
            rolling_beta[i] <- cov_val / var_bench
          }
        }
      }
      rolling_beta[1:59] <- NA

      plot(x$daily$dates, rolling_beta, type = "l",
           col = "orange", lwd = 2,
           main = "Rolling 60-Day Beta",
           xlab = "Date", ylab = "Beta")
      abline(h = x$benchmark$beta, lty = 2)
      grid()
    }
  }

  invisible(x)
}
