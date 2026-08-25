args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else "r_results.csv"
suppressPackageStartupMessages(library(PerformanceAnalytics))
source(file.path("..", "common", "read_asset_prices_binary.R"))
prices <- read_asset_prices_binary(file.path("..", "..", "asset_class_etf_prices.bin"))
rets <- apply(prices$prices, 2L, function(x) diff(x)/head(x, -1L))
spy <- xts::xts(rets[, "SPY"], order.by=prices$dates[-1L])
efa <- xts::xts(rets[, "EFA"], order.by=prices$dates[-1L])
rows <- list()
add <- function(name, fun, reps=100L, atol=1e-10, rtol=1e-8) {
  tm <- system.time(for (i in seq_len(reps)) value <- as.numeric(fun()))[["elapsed"]]
  rows[[length(rows)+1L]] <<- data.frame(case=name, value=value, seconds=tm,
    abs_tol=atol, rel_tol=rtol)
}
add("etf_spy_annualized_return", function() Return.annualized(spy, scale=252))
add("etf_spy_annualized_sd", function() StdDev.annualized(spy, scale=252))
add("etf_spy_sharpe", function() SharpeRatio.annualized(spy, Rf=0, scale=252,
  geometric=FALSE))
add("etf_spy_sortino", function() SortinoRatio(spy, MAR=0))
add("etf_spy_max_drawdown", function() maxDrawdown(spy, invert=TRUE))
add("etf_spy_historical_var95", function() VaR(spy, p=.95, method="historical",
  invert=FALSE))
add("etf_spy_historical_es95", function() ES(spy, p=.95, method="historical",
  invert=FALSE))
add("etf_efa_capm_alpha", function() CAPM.alpha(efa, spy, Rf=0))
add("etf_efa_capm_beta", function() CAPM.beta(efa, spy, Rf=0))
add("etf_efa_capm_beta_bull", function() CAPM.beta.bull(efa, spy, Rf=0))
add("etf_efa_spy_correlation", function() as.numeric(cor(efa,spy)))
add("etf_efa_tracking_error", function() TrackingError(efa, spy, scale=252))
write.csv(do.call(rbind, rows), out, row.names=FALSE, quote=FALSE)
