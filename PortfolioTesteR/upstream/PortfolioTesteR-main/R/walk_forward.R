
#' Walk-Forward Optimization Analysis
#'
#' @description
#' Runs rolling IS/OOS optimization, reselects params each window, and
#' backtests OOS performance (optionally with warmup tails).
#'
#' @param prices Data frame with Date column and symbol columns
#' @param grid   Data frame OR named list; each row/combination is a parameter set
#' @param builder Function(prices, params, ...) -> weights data.frame (Date + assets)
#' @param metric Function(backtest_result) -> scalar score (higher is better).
#'               Defaults to \code{metric_sharpe} if omitted/NULL.
#' @param is_periods Integer, number of in-sample periods
#' @param oos_periods Integer, number of out-of-sample periods
#' @param step Integer, step size for rolling windows (default = oos_periods)
#' @param warmup_periods Integer, warmup periods appended before each OOS
#' @param verbose Logical, print progress
#' @param light_mode Logical, passed to run_param_grid (kept for compatibility)
#' @param precompute_all Logical, precompute indicators once and slice per window
#' @param builder_args List, extra args passed to builder (e.g., indicator_cache)
#' @param n_cores Integer (kept for API compatibility; ignored here)
#'
#' @return An object of class \code{wf_optimization_result}.
#' @export
run_walk_forward <- function(prices,
                             grid,
                             builder,
                             metric = NULL,
                             is_periods = 52,
                             oos_periods = 13,
                             step = NULL,
                             warmup_periods = 0,
                             verbose = FALSE,
                             light_mode = TRUE,
                             precompute_all = TRUE,
                             builder_args = list(),
                             n_cores = 1) {

  # Default metric
  if (is.null(metric)) metric <- metric_sharpe
  if (!is.function(metric)) stop("'metric' must be a function", call. = FALSE)

  # Normalize grid: allow named list OR data.frame
  if (is.list(grid) && !is.data.frame(grid)) {
    grid <- do.call(expand.grid, c(grid, stringsAsFactors = FALSE))
  }
  if (!is.data.frame(grid)) {
    stop("'grid' must be a data frame or a named list convertible via expand.grid()", call. = FALSE)
  }

  # Validation
  if (!is.data.frame(prices)) stop("'prices' must be a data frame", call. = FALSE)
  if (!is.function(builder))  stop("'builder' must be a function", call. = FALSE)
  if (is_periods < 10)        stop("'is_periods' must be at least 10", call. = FALSE)
  if (oos_periods < 2)        stop("'oos_periods' must be at least 2", call. = FALSE)

  # Guard: prevent cadence knobs in grid
  .pt_check_grid_params(grid, "walk-forward optimization")

  # Step default (no overlap)
  if (is.null(step)) step <- oos_periods
  if (step < oos_periods) {
    warning("Step size smaller than OOS period - results will overlap!", immediate. = TRUE)
  }

  # Window splits
  n_rows  <- nrow(prices)
  windows <- .wf_make_splits(n_rows, is_periods, oos_periods, step)
  if (length(windows) == 0) {
    stop(sprintf("Data too short for walk-forward: need at least %d periods, have %d",
                 is_periods + oos_periods, n_rows), call. = FALSE)
  }

  # Optional precompute (e.g., momentum cache)
  if (precompute_all && "lookback" %in% names(grid)) {
    unique_lookbacks <- unique(grid$lookback)
    cache <- list()
    if (exists("calc_momentum", mode = "function", inherits = TRUE)) {
      for (lb in unique_lookbacks) {
        cache[[paste0("momentum_", lb)]] <- calc_momentum(data = prices, lookback = lb)
      }
    }
    if (length(cache) > 0) builder_args$indicator_cache <- cache
  }

  # Header
  if (verbose) {
    freq_factor <- .pt_detect_frequency(prices$Date)
    freq_name <- if (freq_factor == 252) "daily" else if (freq_factor == 52) "weekly" else if (freq_factor == 12) "monthly" else "quarterly"
    cat("=========================================\n")
    cat("       WALK-FORWARD OPTIMIZATION        \n")
    cat("=========================================\n")
    cat("Data:", n_rows, "periods,", freq_name, "frequency\n")
    cat("Windows: IS =", is_periods, "periods, OOS =", oos_periods, "periods\n")
    cat("Step:", step, "periods\n")
    cat("Number of windows:", length(windows), "\n")
    cat("Parameter grid size:", nrow(grid), "\n")
    cat("Mode:", if (light_mode) "Light (fast)" else "Full", "\n\n")
  }

  # Storage
  oos_results <- list()
  chosen_params <- list()
  optimization_details <- list()

  # Loop windows
  for (i in seq_along(windows)) {
    window <- windows[[i]]

    if (verbose) {
      cat(sprintf("Window %d/%d:\n", i, length(windows)))
      is_dates  <- prices$Date[window$is_idx]
      oos_dates <- prices$Date[window$oos_idx]
      cat(sprintf("  IS:  %s to %s (%d periods)\n",
                  format(min(is_dates), "%Y-%m-%d"),
                  format(max(is_dates), "%Y-%m-%d"),
                  length(window$is_idx)))
      cat(sprintf("  OOS: %s to %s (%d periods)\n",
                  format(min(oos_dates), "%Y-%m-%d"),
                  format(max(oos_dates), "%Y-%m-%d"),
                  length(window$oos_idx)))
    }

    # IS slice
    is_prices <- prices[window$is_idx, , drop = FALSE]

    # IS cache slice
    is_builder_args <- builder_args
    if (precompute_all && !is.null(builder_args$indicator_cache)) {
      is_cache <- list()
      for (nm in names(builder_args$indicator_cache)) {
        full_ind <- builder_args$indicator_cache[[nm]]
        is_cache[[nm]] <- full_ind[window$is_idx, , drop = FALSE]
      }
      is_builder_args$indicator_cache <- is_cache
    }

    # Optimize on IS
    if (verbose) cat("  Optimizing on in-sample data...\n")
    t0 <- proc.time()
    opt_result <- run_param_grid(
      prices = is_prices,
      grid = grid,
      builder = builder,
      metric = metric,
      name_prefix = paste0("IS_", i),
      verbose = FALSE,
      light_mode = light_mode,
      precompute_returns = TRUE,
      builder_args = is_builder_args,
      n_cores = n_cores
    )

    if (is.null(opt_result$best_params)) {
      if (verbose) cat("  WARNING: No valid IS results, skipping window\n\n")
      next
    }

    chosen_params[[i]] <- opt_result$best_params

    if (verbose) {
      cat(sprintf("  Best IS parameters: %s\n",
                  paste(sprintf("%s=%s",
                                names(opt_result$best_params),
                                opt_result$best_params),
                        collapse = ", ")))
      cat(sprintf("  Best IS score: %.4f\n", opt_result$optimization_score))
    }

    # Build OOS (with/without warmup) + VALIDATE weights against the
    # exact prices slice used for the forthcoming backtest
    if (warmup_periods > 0 && warmup_periods < length(window$is_idx)) {
      # Extended window with warmup tail from IS
      warmup_idx   <- tail(window$is_idx, warmup_periods)
      extended_idx <- c(warmup_idx, window$oos_idx)
      extended_pr  <- prices[extended_idx, , drop = FALSE]

      # OOS cache slice for extended window
      oos_builder_args <- builder_args
      if (precompute_all && !is.null(builder_args$indicator_cache)) {
        oos_cache <- list()
        for (nm in names(builder_args$indicator_cache)) {
          full_ind <- builder_args$indicator_cache[[nm]]
          oos_cache[[nm]] <- full_ind[extended_idx, , drop = FALSE]
        }
        oos_builder_args$indicator_cache <- oos_cache
      }

      # Build weights on extended window
      extended_w <- if (length(oos_builder_args) > 0) {
        do.call(builder, c(list(extended_pr, opt_result$best_params), oos_builder_args))
      } else {
        builder(extended_pr, opt_result$best_params)
      }

      # [OK] Validate against extended window
      extended_w <- .pt_validate_weights(
        extended_pr, extended_w,
        context = sprintf("WF window %d (extended)", i),
        allow_sparse = TRUE
      )

      # Drop warmup rows, then validate against pure OOS
      oos_start   <- length(warmup_idx) + 1
      oos_weights <- extended_w[oos_start:nrow(extended_w), , drop = FALSE]
      oos_prices  <- prices[window$oos_idx, , drop = FALSE]

      # [OK] Validate OOS alignment exactly
      oos_weights <- .pt_validate_weights(
        oos_prices, oos_weights,
        context = sprintf("WF window %d (OOS from warmup)", i),
        allow_sparse = TRUE
      )

    } else {
      # No warmup
      oos_prices <- prices[window$oos_idx, , drop = FALSE]

      # OOS cache slice for OOS window
      oos_builder_args <- builder_args
      if (precompute_all && !is.null(builder_args$indicator_cache)) {
        oos_cache <- list()
        for (nm in names(builder_args$indicator_cache)) {
          full_ind <- builder_args$indicator_cache[[nm]]
          oos_cache[[nm]] <- full_ind[window$oos_idx, , drop = FALSE]
        }
        oos_builder_args$indicator_cache <- oos_cache
      }

      # Build weights on OOS window
      oos_weights <- if (length(oos_builder_args) > 0) {
        do.call(builder, c(list(oos_prices, opt_result$best_params), oos_builder_args))
      } else {
        builder(oos_prices, opt_result$best_params)
      }

      # [OK] Validate OOS weights vs OOS prices
      oos_weights <- .pt_validate_weights(
        oos_prices, oos_weights,
        context = sprintf("WF window %d (OOS)", i),
        allow_sparse = TRUE
      )
    }

    # OOS backtest (full mode for stitching)
    oos_bt <- run_backtest(
      prices = oos_prices,
      weights = oos_weights,
      initial_capital = 100000,
      name = paste0("OOS_", i),
      verbose = FALSE
    )

    oos_score <- metric(oos_bt)
    elapsed   <- unname((proc.time() - t0)["elapsed"])

    if (verbose) {
      cat(sprintf("  OOS score: %.4f\n", oos_score))
      cat(sprintf("  OOS return: %.1f%%\n", oos_bt$total_return * 100))
      cat(sprintf("  Window time: %.1fs\n\n", elapsed))
    }

    # Store results
    oos_results[[i]] <- oos_bt
    optimization_details[[i]] <- list(
      window = i,
      is_start = min(prices$Date[window$is_idx]),
      is_end   = max(prices$Date[window$is_idx]),
      oos_start = min(prices$Date[window$oos_idx]),
      oos_end   = max(prices$Date[window$oos_idx]),
      is_score  = opt_result$optimization_score,
      oos_score = oos_score,
      oos_return = oos_bt$total_return,
      n_transactions = if (!is.null(oos_bt$n_transactions)) oos_bt$n_transactions else NA,
      optimization_time = elapsed,
      params = opt_result$best_params
    )
  }

  # Stitch & summarize
  if (verbose) cat("Stitching out-of-sample results...\n")
  oos_stitched <- wf_stitch(oos_results)

  summary_df <- do.call(
    rbind.data.frame,
    lapply(optimization_details, function(x) {
      if (is.null(x)) return(NULL)
      params_df <- as.data.frame(x$params)
      data.frame(
        Window     = x$window,
        IS_Score   = x$is_score,
        OOS_Score  = x$oos_score,
        OOS_Return = x$oos_return,
        params_df,
        stringsAsFactors = FALSE
      )
    })
  )

  structure(
    list(
      oos_stitched = oos_stitched,
      oos_results  = oos_results,
      chosen_params = chosen_params,
      windows = windows,
      optimization_summary  = summary_df,
      optimization_details  = optimization_details,
      grid = grid,
      is_periods  = is_periods,
      oos_periods = oos_periods,
      step = step
    ),
    class = "wf_optimization_result"
  )
}

#' Create Window Splits for Walk-Forward
#' @keywords internal
.wf_make_splits <- function(n, is_periods, oos_periods, step) {
  windows <- list()
  start_idx <- 1
  k <- 1
  while (start_idx + is_periods + oos_periods - 1 <= n) {
    is_end    <- start_idx + is_periods - 1
    oos_start <- is_end + 1
    oos_end   <- oos_start + oos_periods - 1
    windows[[k]] <- list(is_idx = start_idx:is_end, oos_idx = oos_start:oos_end)
    k <- k + 1
    start_idx <- start_idx + step
  }
  windows
}

#' Stitch Out-of-Sample Results (overlap-safe)
#'
#' @description
#' Concatenates OOS backtests and safely compounds returns on overlapping dates.
#'
#' @param oos_results List of backtest_result objects, each with $portfolio_values and $dates.
#' @param initial_value Numeric starting value for the stitched equity curve (default 100000).
#'
#' @return Data frame with columns: \code{Date}, \code{Value}.
#' @importFrom stats aggregate
#' @export
wf_stitch <- function(oos_results, initial_value = 100000) {
  oos_results <- Filter(Negate(is.null), oos_results)
  if (!length(oos_results)) {
    return(data.frame(Date = as.Date(character(0)), Value = numeric(0)))
  }

  all_returns <- numeric()
  all_dates   <- as.Date(character())

  for (bt in oos_results) {
    pv <- bt$portfolio_values
    if (length(pv) > 1) {
      period_returns <- diff(pv) / head(pv, -1)
      period_dates   <- as.Date(bt$dates[-1])
      all_returns <- c(all_returns, period_returns)
      all_dates   <- c(all_dates, period_dates)
    }
  }

  if (!length(all_returns)) {
    return(data.frame(Date = as.Date(character(0)), Value = numeric(0)))
  }

  # If overlap was allowed, compound same-date returns
  if (any(duplicated(all_dates))) {
    df <- data.frame(Date = all_dates, Return = all_returns)
    # compound returns per date: prod(1+r_i) - 1
    df <- aggregate(Return ~ Date, df, FUN = function(x) prod(1 + x) - 1)
    all_dates   <- df$Date
    all_returns <- df$Return
  }

  # Sort and cumulate
  ord <- order(all_dates)
  all_dates   <- all_dates[ord]
  all_returns <- all_returns[ord]

  values <- initial_value * cumprod(1 + all_returns)
  data.frame(
    Date  = c(min(all_dates), all_dates),
    Value = c(initial_value, values),
    stringsAsFactors = FALSE
  )
}

#' Generate Walk-Forward Report
#'
#' @description
#' Prints a concise summary of a \code{wf_optimization_result}: configuration,
#' stitched OOS performance, and parameter stability.
#'
#' @param wf A \code{wf_optimization_result} object (from \code{run_walk_forward()}).
#' @param digits Integer; number of digits when printing numeric values (default 4).
#'
#' @return Invisibly returns the optimization summary data frame.
#' @export
wf_report <- function(wf, digits = 4) {
  cat("=========================================\n")
  cat("   WALK-FORWARD OPTIMIZATION REPORT     \n")
  cat("=========================================\n\n")
  cat("Configuration:\n")
  cat("--------------\n")
  cat("Total windows:", length(wf$windows), "\n")
  cat("IS periods:", wf$is_periods, "\n")
  cat("OOS periods:", wf$oos_periods, "\n")
  cat("Step size:", .pt_null_coalesce(wf$step, wf$oos_periods), "\n")
  cat("Parameter combinations tested per window:", nrow(wf$grid), "\n\n")

  if (!is.null(wf$oos_stitched) && nrow(wf$oos_stitched) > 1) {
    total_return <- (tail(wf$oos_stitched$Value, 1) / head(wf$oos_stitched$Value, 1) - 1)
    cat("Overall Performance:\n")
    cat("-------------------\n")
    cat(sprintf("Total stitched return: %.1f%%\n", total_return * 100))

    date_range <- as.numeric(difftime(
      max(wf$oos_stitched$Date),
      min(wf$oos_stitched$Date),
      units = "days"
    ))
    years <- date_range / 365.25
    if (years > 0) {
      cagr <- (1 + total_return)^(1 / years) - 1
      cat(sprintf("Annualized return: %.1f%%\n", cagr * 100))
      cat(sprintf("Period: %.1f years\n", years))
    }
  }

  if (length(wf$chosen_params) > 1) {
    cat("\nParameter Stability:\n")
    cat("-------------------\n")
    all_params <- unique(unlist(lapply(wf$chosen_params, names)))
    for (param in all_params) {
      values <- sapply(wf$chosen_params, function(x) x[[param]])
      unique_vals <- unique(values[!is.na(values)])
      if (length(unique_vals)) {
        tab <- table(values)
        most_common <- names(tab)[which.max(tab)]
        cat(sprintf("%s: Most common = %s, Unique values = %d\n",
                    param, most_common, length(unique_vals)))
      }
    }
  }

  invisible(wf$optimization_summary)
}

#' Print a wf_optimization_result
#'
#' @param x A `wf_optimization_result` object returned by [run_walk_forward()].
#' @param ... Additional arguments passed to methods (ignored).
#'
#' @return Invisibly returns `x`.
#' @export
#' @method print wf_optimization_result
#' @aliases print.wf_optimization_result
print.wf_optimization_result <- function(x, ...) {
  cat("Walk-Forward Optimization Result\n")
  cat("=================================\n")
  cat("Windows:", length(x$windows), "\n")
  cat("IS/OOS/Step:", x$is_periods, "/", x$oos_periods, "/", .pt_null_coalesce(x$step, x$oos_periods), "\n")
  if (!is.null(x$oos_stitched) && nrow(x$oos_stitched) > 1) {
    total_return <- (tail(x$oos_stitched$Value, 1) / head(x$oos_stitched$Value, 1) - 1)
    cat(sprintf("Total OOS return: %.1f%%\n", total_return * 100))
  }
  cat("\nUse wf_report() for detailed results\n")
  invisible(x)
}


#' Plot Walk-Forward Results
#'
#' Visual diagnostics for a \code{wf_optimization_result} returned by
#' [run_walk_forward()]. Supported types:
#' \itemize{
#'   \item "parameters": best parameter values chosen per window.
#'   \item "is_oos": in-sample vs out-of-sample scores by window.
#'   \item "equity": stitched out-of-sample equity curve.
#'   \item "drawdown": drawdown of the stitched OOS curve.
#'   \item "windows": per-window bar/line chart of an OOS metric (see \code{metric}).
#'   \item "stability": summary of parameter stability.
#'   \item "distributions": distributions of IS/OOS scores across windows.
#' }
#'
#' @param x A \code{wf_optimization_result} from [run_walk_forward()].
#' @param y Ignored.
#' @param type One of "parameters","is_oos","equity","drawdown",
#'   "windows","stability","distributions".
#' @param param Character vector of parameter names to include for
#'   "parameters"/"stability"/"distributions". If \code{NULL}, uses all.
#' @param metric Character; column to plot for "is_oos" or "windows"
#'   (e.g., "OOS_Return" or "OOS_Score"). Ignored for other types.
#' @param main,sub,xlab,ylab Base plotting annotations.
#' @param ... Additional plot options (type-specific).
#'
#' @return Invisibly, a small list describing the plot.
#'
#' @examples
#' \donttest{
#'   data(sample_prices_weekly)
#'   b <- function(prices, params, ...) {
#'     weight_equally(filter_top_n(calc_momentum(prices, params$lookback),
#'                                 params$top_n))
#'   }
#'   wf <- run_walk_forward(
#'     prices = sample_prices_weekly,
#'     grid   = list(lookback = c(8, 12, 26), top_n = c(5, 10)),
#'     builder = b,
#'     is_periods = 52, oos_periods = 26, step = 26
#'   )
#'   plot(wf, type = "parameters")
#'   plot(wf, type = "is_oos", metric = "OOS_Score")
#' }
#'
#' @seealso [run_walk_forward()], [wf_report()], [print.wf_optimization_result()]
#'
#' @method plot wf_optimization_result
#' @export
#' @importFrom graphics axis box layout plot.new points title
plot.wf_optimization_result <- function(
    x, y = NULL,
    type = c("parameters","is_oos","equity","drawdown","windows","stability","distributions"),
    param = NULL,          # which parameter(s) for parameters/stability/distributions
    metric = NULL,         # for is_oos OR windows (choose which cols)
    main = NULL, sub = NULL, xlab = NULL, ylab = NULL, ...
) {
  type <- match.arg(type)

  # par safety
  if (!exists(".with_par", mode = "function")) {
    .with_par <- function(expr) { old <- par(no.readonly = TRUE); on.exit(par(old)); force(expr) }
  }

  os <- x$optimization_summary
  if (is.null(os) || !nrow(os)) stop("wf_optimization_result has empty optimization_summary.")

  # helper to detect parameter columns from summary
  fixed_cols <- c("Window","IS_Score","OOS_Score","OOS_Return",
                  "n_transactions","optimization_time")
  param_cols <- setdiff(names(os), intersect(names(os), fixed_cols))

  # ---------------- PARAMETERS (evolution across windows) ----------------
  # ---- PARAMETERS (evolution across windows) ----
  if (type == "parameters") {
    # always use a simple 1..n sequence for x (avoids list/Date issues)
    win <- seq_len(nrow(os))

    # fields that aren't parameters
    fixed_cols <- c("Window","IS_Score","OOS_Score","OOS_Return",
                    "n_transactions","optimization_time")
    param_cols <- setdiff(names(os), fixed_cols)
    if (!is.null(param)) param_cols <- intersect(param_cols, param)
    if (!length(param_cols)) stop("No parameter columns to plot.")

    k <- length(param_cols); r <- ceiling(sqrt(k)); c <- ceiling(k / r)

    .with_par({
      par(mfrow = c(r, c), mar = c(4, 4, 3, 1) + 0.1)
      for (p in param_cols) {
        v <- os[[p]]
        if (is.list(v)) v <- unlist(v, use.names = FALSE)     # flatten list-cols

        vn <- suppressWarnings(as.numeric(v))
        if (all(is.finite(vn))) {
          plot(win, vn, type = "b", pch = 19,
               xlab = if (is.null(xlab)) "Window" else xlab,
               ylab = if (is.null(ylab)) p        else ylab,
               main = if (is.null(main)) sprintf("Parameter: %s", p) else main, ...)
          grid()
        } else {
          f  <- factor(v, levels = unique(v))
          yy <- as.integer(f)
          plot(win, yy, type = "b", pch = 19, yaxt = "n",
               xlab = if (is.null(xlab)) "Window" else xlab,
               ylab = if (is.null(ylab)) p        else ylab,
               main = if (is.null(main)) sprintf("Parameter: %s", p) else main, ...)
          axis(2, at = seq_along(levels(f)), labels = levels(f), las = 2)
          grid()
        }
      }
    })

    return(invisible(list(kind="parameters", params=param_cols, n=k)))
  }


  # ---------------- IS vs OOS SCATTER ----------------
  if (type == "is_oos") {
    is_col  <- if (!is.null(metric) && length(metric) >= 1) metric[1] else
      if ("IS_Score"  %in% names(os))  "IS_Score"  else names(os)[2]
    oos_col <- if (!is.null(metric) && length(metric) >= 2) metric[2] else
      if ("OOS_Score" %in% names(os)) "OOS_Score" else names(os)[3]
    stopifnot(all(c(is_col, oos_col) %in% names(os)))
    IS <- os[[is_col]]; OOS <- os[[oos_col]]
    ok <- is.finite(IS) & is.finite(OOS); if (!any(ok)) stop("No finite IS/OOS pairs.")
    rng <- range(c(IS[ok], OOS[ok]), na.rm = TRUE)
    rho <- suppressWarnings(cor(IS, OOS, use = "complete.obs"))
    .with_par({
      plot(IS, OOS, xlab = is_col, ylab = oos_col,
           main = if (is.null(main)) "IS vs OOS" else main,
           xlim = rng, ylim = rng, pch = 19, ...); abline(0,1,lty=2,col="gray40"); grid()
      legend("topleft", legend = sprintf("rho = %.2f", rho), bty = "n")
    })
    return(invisible(list(kind="is_oos", is_col=is_col, oos_col=oos_col, cor=rho)))
  }

  # ---------------- EQUITY (stitched OOS curve) ----------------
  if (type == "equity") {
    eq <- x$oos_stitched
    if (is.null(eq) || !nrow(eq)) stop("oos_stitched is empty.")
    .with_par({
      plot(eq$Date, eq$Value, type="l", lwd=2,
           xlab = if (is.null(xlab)) "Date" else xlab,
           ylab = if (is.null(ylab)) "Portfolio Value" else ylab,
           main = if (is.null(main)) "Walk-Forward Equity Curve" else main, ...)
      grid()
    })
    return(invisible(list(kind="equity", n=nrow(eq))))
  }

  # ---------------- DRAWDOWN (from stitched OOS) ----------------
  if (type == "drawdown") {
    eq <- x$oos_stitched
    if (is.null(eq) || !nrow(eq)) stop("oos_stitched is empty.")
    v  <- eq$Value; dd <- 1 - v / cummax(v)
    .with_par({
      plot(eq$Date, dd * 100, type="l", lwd=2,
           xlab = if (is.null(xlab)) "Date" else xlab,
           ylab = if (is.null(ylab)) "Drawdown (%)" else ylab,
           main = if (is.null(main)) "Walk-Forward Drawdown" else main, ...)
      abline(h=0,lty=2,col="gray40"); grid()
    })
    return(invisible(list(kind="drawdown", min_dd = min(dd, na.rm=TRUE))))
  }

  # ---------------- WINDOWS (per-window bars) ----------------
  if (type == "windows") {
    win <- if ("Window" %in% names(os)) os$Window else seq_len(nrow(os))
    # choose metric for bars: default OOS_Return if present, else OOS_Score
    mcol <- if (!is.null(metric)) metric else if ("OOS_Return" %in% names(os)) "OOS_Return" else "OOS_Score"
    stopifnot(mcol %in% names(os))
    vals <- os[[mcol]]
    .with_par({
      par(mar = c(5, 4, 4, 2) + 0.1)
      bp <- barplot(height = vals, names.arg = win,
                    xlab = if (is.null(xlab)) "Window" else xlab,
                    ylab = if (is.null(ylab)) mcol else ylab,
                    main = if (is.null(main)) sprintf("Per-window %s", mcol) else main, ...)
      abline(h = 0, col = "gray40", lty = 2); grid(nx = 0, ny = NULL)
    })
    return(invisible(list(kind="windows", metric=mcol)))
  }

  # ---------------- STABILITY (freq; 1D or 2D) ----------------
  if (type == "stability") {
    if (is.null(param) || length(param) == 0) {
      # pick first two params if not specified
      if (!length(param_cols)) stop("No parameter columns in optimization_summary.")
      param <- head(param_cols, 2)
    }
    if (length(param) == 1) {
      v <- os[[param]]; if (is.list(v)) v <- unlist(v, use.names = FALSE)
      f <- factor(v)
      .with_par({
        barplot(table(f), main = if (is.null(main)) sprintf("Stability: %s", param) else main,
                xlab = param, ylab = "Frequency", ...)
        grid(nx = 0, ny = NULL)
      })
      return(invisible(list(kind="stability", params=param, dim=1)))
    } else {
      # ... inside plot.wf_optimization_result(), in the stability 2D else-branch:

      p1 <- param[1]; p2 <- param[2]
      a <- os[[p1]]; b <- os[[p2]]
      if (is.list(a)) a <- unlist(a, use.names = FALSE)
      if (is.list(b)) b <- unlist(b, use.names = FALSE)

      A <- factor(a); B <- factor(b)
      tab <- table(B, A)                     # rows = B (y), cols = A (x)
      z   <- as.matrix(tab)

      xlabs <- colnames(z);  ylabs <- rownames(z)
      xvals <- seq_len(ncol(z)); yvals <- seq_len(nrow(z))

      pal <- grDevices::hcl.colors(100, "YlOrRd", rev = TRUE)

      .with_par({
        # transpose so nrow(zT)=length(x), ncol(zT)=length(y)
        zT <- t(z)
        image(xvals, yvals, zT, col = pal,
              xlab = p1, ylab = p2,
              main = if (is.null(main)) sprintf("Stability: %s x %s", p1, p2) else main, ...)
        axis(1, at = xvals, labels = xlabs, las = 2)
        axis(2, at = yvals, labels = ylabs, las = 2)
        box()
      })
      return(invisible(list(kind="stability", params=c(p1,p2), dim=2)))

    }
  }

  # ---------------- DISTRIBUTIONS (per-parameter) ----------------
  if (type == "distributions") {
    if (!length(param_cols)) stop("No parameter columns in optimization_summary.")
    if (!is.null(param)) param_cols <- intersect(param_cols, param)
    k <- length(param_cols); r <- ceiling(sqrt(k)); c <- ceiling(k / r)
    .with_par({
      par(mfrow = c(r, c), mar = c(4, 4, 3, 1) + 0.1)
      for (p in param_cols) {
        v <- os[[p]]; if (is.list(v)) v <- unlist(v, use.names = FALSE)
        vn <- suppressWarnings(as.numeric(v))
        if (all(is.finite(vn))) {
          hist(vn, breaks = "FD",
               main = if (is.null(main)) sprintf("Distribution: %s", p) else main,
               xlab = p, ylab = "Count", col = "gray85", border = "white", ...)
          grid()
        } else {
          barplot(table(v),
                  main = if (is.null(main)) sprintf("Distribution: %s", p) else main,
                  xlab = p, ylab = "Count", ...)
          grid(nx = 0, ny = NULL)
        }
      }
    })
    return(invisible(list(kind="distributions", params=param_cols)))
  }

  message(sprintf("Plot type '%s' not implemented.", type))
  invisible(list(type = type))
}

