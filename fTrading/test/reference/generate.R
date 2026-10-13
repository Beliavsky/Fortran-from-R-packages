# SPDX-License-Identifier: LGPL-2.0-or-later
# Deterministic oracle fixtures for the Fortran translation of fTrading 3042.79.
# Run from the fTrading package root: Rscript test/reference/generate.R
# Only the genuine installed package supplies expected numerical results.
# No functions are patched or reimplemented. Original timeSeries paths are used
# where the upstream numeric-vector interface raises an error.

suppressPackageStartupMessages(library(fTrading))
stopifnot(as.character(packageVersion("fTrading")) == "3042.79")
ft <- function(name) getExportedValue("fTrading", name)

# Confirm the installed oracle matches the bundled, unmodified upstream source.
source_environment <- new.env(parent = asNamespace("fTrading"))
for (file in c("BenchmarkAnalysis.R", "TechnicalAnalysis.R", "rollFun.R")) {
    sys.source(file.path("upstream", "R", file), envir = source_environment)
}
original_functions <- ls(source_environment, all.names = TRUE)
for (name in original_functions) {
    original <- source_environment[[name]]
    installed <- ft(name)
    stopifnot(identical(formals(original), formals(installed)))
    stopifnot(identical(deparse(body(original)), deparse(body(installed))))
}
stopifnot(length(original_functions) == 45L)

n <- 96L
i <- seq_len(n)
close <- 100 + 0.09 * i + 3 * sin(0.47 * i) + 1.4 * cos(0.19 * i)
close[c(14, 33)] <- close[c(13, 32)]
open <- close - 0.4 + 0.2 * cos(0.23 * i)
high <- pmax(open, close) + 1.2 + 0.35 * (i %% 4)
low <- pmin(open, close) - 0.9 - 0.2 * (i %% 5)
volume <- 1000 + 25 * (i %% 9) + 200 * sin(0.51 * i)
volume[20] <- volume[19]
levels <- log(close)
rates <- 0.0001 * ((seq_len(n - 1L) %% 5) - 2)
missing <- close[1:16]
missing[c(3, 11)] <- NA_real_
positions <- c(-1, -1, 1, 1, 1, -1, -1, 1, -1, -1, 1, 1)
returns <- c(-0.03, 0.02, 0.04, -0.01, 0, 0.02, -0.05, 0.03, -0.02, 0.04, 0.01, 0)
positions_missing <- positions
positions_missing[4] <- NA_real_
returns_missing <- returns
returns_missing[8] <- NA_real_
dates <- as.character(as.Date("2020-01-01") + seq_len(n) - 1L)
ts_high <- timeSeries(high, charvec = dates, units = "High")
ts_low <- timeSeries(low, charvec = dates, units = "Low")
ts_close <- timeSeries(close, charvec = dates, units = "Close")
ts_ohlc <- timeSeries(cbind(Open = open, High = high, Low = low, Close = close), charvec = dates)
ts_positions <- timeSeries(positions, charvec = dates[seq_along(positions)], units = "Position")

# Explicitly verify that the original vector-path defects are still present.
vector_failures <- c(
    chaikinvTA = tryCatch({ft("chaikinvTA")(high, low); "unexpected success"}, error = conditionMessage),
    williamsadTA = tryCatch({ft("williamsadTA")(high, low, close); "unexpected success"}, error = conditionMessage),
    williamsrTA = tryCatch({ft("williamsrTA")(high, low, close); "unexpected success"}, error = conditionMessage)
)
stopifnot(all(vector_failures != "unexpected success"))

fixtures <- list()
add <- function(name, value) {
    stopifnot(is.null(fixtures[[name]]))
    fixtures[[name]] <<- as.vector(value)
}
for (name in c("open", "high", "low", "close", "volume", "levels", "rates", "missing",
               "positions", "returns", "positions_missing", "returns_missing")) {
    add(paste0("input_", name), get(name))
}
add("sharpeRatio_default", ft("sharpeRatio")(levels))
add("sharpeRatio_scalar", ft("sharpeRatio")(levels, r = 0.003, scale = 7))
add("sharpeRatio_rates", ft("sharpeRatio")(levels, r = rates, scale = 3))
add("sharpeRatio_short", ft("sharpeRatio")(levels[1:2]))
add("sterlingRatio_default", ft("sterlingRatio")(levels))
add("sterlingRatio_rising", ft("sterlingRatio")(c(1, 2, 4, 8)))
drawdown <- ft("maxDrawDown")(levels)
add("maxDrawDown_levels_value", drawdown$maxdrawdown)
add("maxDrawDown_levels_from", drawdown$from)
add("maxDrawDown_levels_to", drawdown$to)
drawdown <- ft("maxDrawDown")(c(10, 12, 9, 12, 9, 14, 11, 14, 11))
add("maxDrawDown_ties_value", drawdown$maxdrawdown)
add("maxDrawDown_ties_from", drawdown$from)
add("maxDrawDown_ties_to", drawdown$to)
drawdown <- ft("maxDrawDown")(c(1, 2, 4, 8))
add("maxDrawDown_rising_value", drawdown$maxdrawdown)
add("maxDrawDown_rising_from", drawdown$from)
add("maxDrawDown_rising_to", drawdown$to)
add("rollFun_mean", ft("rollFun")(close, n = 7, FUN = mean))
add("rollFun_sum", ft("rollFun")(close, n = 7, FUN = sum))
add("rollFun_min", ft("rollFun")(close, n = 7, FUN = min))
add("rollFun_max", ft("rollFun")(close, n = 7, FUN = max))
add("rollFun_var", ft("rollFun")(close, n = 7, FUN = var))
add("rollFun_sd", ft("rollFun")(close, n = 7, FUN = sd))
add("rollFun_padded", ft("rollFun")(close, n = 4, trim = FALSE, FUN = mean))
add("rollFun_na_windows", ft("rollFun")(missing, n = 3, trim = FALSE, FUN = mean))
add("rollFun_na_removed", ft("rollFun")(missing, n = 3, trim = FALSE, na.rm = TRUE, FUN = mean))
add("rollFun_callback", ft("rollFun")(close, n = 6, FUN = function(window) sum(window^2)))
add("rollFun_unit", ft("rollFun")(close, n = 1, FUN = sum))
add("rollVar_default", ft("rollVar")(close))
add("rollVar_biased", ft("rollVar")(close, n = 7, trim = FALSE, unbiased = FALSE))
add("rollVar_na_removed", ft("rollVar")(missing, n = 3, na.rm = TRUE, unbiased = FALSE))
add("rollVar_unit", ft("rollVar")(close, n = 1))
add("SMA_default", ft("SMA")(close))
add("SMA_custom", ft("SMA")(close, n = 13))
add("emaTA_default", ft("emaTA")(close))
add("emaTA_zero", ft("emaTA")(close, lambda = 0))
add("emaTA_period", ft("emaTA")(close, lambda = 7))
add("emaTA_one", ft("emaTA")(close, lambda = 1))
add("emaTA_decay", ft("emaTA")(close, lambda = 0.35))
add("emaTA_startup", ft("emaTA")(close, lambda = 0.3, startup = 4))
add("emaTA_oversized", ft("emaTA")(close, lambda = 0.3, startup = 97))
add("emaTA_short", ft("emaTA")(close[1:6]))
add("EWMA_default_startup", ft("EWMA")(close, lambda = 0.4))
add("EWMA_custom_startup", ft("EWMA")(close, lambda = 11, startup = 6))
add("biasTA_default", ft("biasTA")(close))
add("biasTA_custom", ft("biasTA")(close, lag = 3))
add("rocTA_default", ft("rocTA")(close))
add("rocTA_custom", ft("rocTA")(close, lag = 3))
add("momTA_default", ft("momTA")(close))
add("momTA_custom", ft("momTA")(close, lag = 3))
add("medpriceTA_default", ft("medpriceTA")(high, low))
add("vohlTA_default", ft("vohlTA")(high, low))
add("vorTA_default", ft("vorTA")(high, low))
add("typicalpriceTA_default", ft("typicalpriceTA")(high, low, close))
add("wcloseTA_default", ft("wcloseTA")(high, low, close))
add("oscTA_default", ft("oscTA")(close))
add("oscTA_custom", ft("oscTA")(close, lag1 = 4, lag2 = 9))
add("macdTA_default", ft("macdTA")(close))
add("macdTA_custom", ft("macdTA")(close, lag1 = 4, lag2 = 9))
add("cdsTA_default", ft("cdsTA")(close))
add("cdsTA_custom", ft("cdsTA")(close, lag1 = 4, lag2 = 9, lag3 = 3))
add("cdoTA_default", ft("cdoTA")(close))
add("cdoTA_custom", ft("cdoTA")(close, lag1 = 4, lag2 = 9, lag3 = 3))
add("fpkTA_default", ft("fpkTA")(close, high, low))
add("fpkTA_custom", ft("fpkTA")(close, high, low, lag = 3))
add("wprTA_default", ft("wprTA")(close, high, low))
add("wprTA_custom", ft("wprTA")(close, high, low, lag = 3))
add("fpdTA_default", ft("fpdTA")(close, high, low))
add("fpdTA_custom", ft("fpdTA")(close, high, low, lag1 = 3, lag2 = 4))
add("spdTA_default", ft("spdTA")(close, high, low))
add("spdTA_custom", ft("spdTA")(close, high, low, lag1 = 3, lag2 = 4, lag3 = 6))
add("apdTA_default", ft("apdTA")(close, high, low))
add("apdTA_custom", ft("apdTA")(close, high, low, lag1 = 3, lag2 = 4, lag3 = 6, lag4 = 7))
add("stochasticTA_fast", ft("stochasticTA")(close, high, low))
add("stochasticTA_slow", ft("stochasticTA")(close, high, low, lag1 = 3, lag2 = 4, lag3 = 6, type = "slow"))
add("rsiTA_default", ft("rsiTA")(close))
add("rsiTA_custom", ft("rsiTA")(close, lag = 4))
add("accelTA_default", ft("accelTA")(close))
add("accelTA_custom", ft("accelTA")(close, n = 4))
add("adiTA_default", ft("adiTA")(high, low, close, volume))
add("adoscillatorTA_default", ft("adoscillatorTA")(open, high, low, close))
add("bollingerTA_default", ft("bollingerTA")(close))
add("bollingerTA_custom", ft("bollingerTA")(close, lag = 7, n.sd = 1.5))
add("chaikinoTA_default", ft("chaikinoTA")(high, low, close, volume))
add("chaikinoTA_custom", ft("chaikinoTA")(high, low, close, volume, lag1 = 6, lag2 = 4))
add("chaikinvTA_default_ts", ft("chaikinvTA")(ts_high, ts_low))
add("chaikinvTA_custom_ts", ft("chaikinvTA")(ts_high, ts_low, lag1 = 6, lag2 = 4))
add("garmanklassTA_default", ft("garmanklassTA")(open, high, low, close))
add("nviTA_default", ft("nviTA")(close, volume))
add("obvTA_default", ft("obvTA")(close, volume))
add("pviTA_default", ft("pviTA")(close, volume))
add("pvtrendTA_default", ft("pvtrendTA")(close, volume))
add("williamsadTA_ts", ft("williamsadTA")(ts_high, ts_low, ts_close))
add("williamsrTA_default_ts", ft("williamsrTA")(ts_high, ts_low, ts_close))
add("williamsrTA_custom_ts", ft("williamsrTA")(ts_high, ts_low, ts_close, lag = 7))
add("dailyTA_ema", ft(".dailyTA")(ts_ohlc, indicator = "ema", lag = c(6)))
add("dailyTA_bias", ft(".dailyTA")(ts_ohlc, indicator = "bias", lag = c(6)))
add("dailyTA_medprice", ft(".dailyTA")(ts_ohlc, indicator = "medprice", lag = c(6)))
add("dailyTA_typicalprice", ft(".dailyTA")(ts_ohlc, indicator = "typicalprice", lag = c(6)))
add("dailyTA_wclose", ft(".dailyTA")(ts_ohlc, indicator = "wclose", lag = c(6)))
add("dailyTA_roc", ft(".dailyTA")(ts_ohlc, indicator = "roc", lag = c(6)))
add("dailyTA_mom", ft(".dailyTA")(ts_ohlc, indicator = "mom", lag = c(6)))
add("dailyTA_osc", ft(".dailyTA")(ts_ohlc, indicator = "osc", lag = c(4, 9)))
add("dailyTA_macd", ft(".dailyTA")(ts_ohlc, indicator = "macd", lag = c(4, 9)))
add("dailyTA_cds", ft(".dailyTA")(ts_ohlc, indicator = "cds", lag = c(4, 9, 3)))
add("dailyTA_cdo", ft(".dailyTA")(ts_ohlc, indicator = "cdo", lag = c(4, 9, 3)))
add("dailyTA_vohl", ft(".dailyTA")(ts_ohlc, indicator = "vohl", lag = c(6)))
add("dailyTA_vor", ft(".dailyTA")(ts_ohlc, indicator = "vor", lag = c(6)))
add("dailyTA_rsi", ft(".dailyTA")(ts_ohlc, indicator = "rsi", lag = c(6)))
add("dailyTA_fpk", ft(".dailyTA")(ts_ohlc, indicator = "fpk", lag = c(6)))
add("dailyTA_fpd", ft(".dailyTA")(ts_ohlc, indicator = "fpd", lag = c(4, 3)))
add("dailyTA_spd", ft(".dailyTA")(ts_ohlc, indicator = "spd", lag = c(4, 3, 6)))
add("dailyTA_apd", ft(".dailyTA")(ts_ohlc, indicator = "apd", lag = c(4, 3, 6, 7)))
add("dailyTA_wpr", ft(".dailyTA")(ts_ohlc, indicator = "wpr", lag = c(6)))
add("dailyTA_select", ft(".dailyTA")(ts_ohlc, indicator = "ema", select = "Open", lag = 6))
add("hitRate_default", ft(".hitRate")(returns, positions))
add("hitRate_missing", ft(".hitRate")(returns_missing, positions_missing))
signals <- ft(".tradeSignals")(ts_positions)
add("tradeSignals_values", signals)
add("tradeSignals_indices", match(as.numeric(time(signals)), as.numeric(time(ts_positions))))
add("input_signal_times", as.numeric(time(signals)))
add("tradeLengths", ft(".tradeLengths")(signals))

# Simple portable text format: record count, then label, length, scalar rows.
connection <- file("test/reference/fixtures.txt", open = "wt")
writeLines(as.character(length(fixtures)), connection)
for (name in names(fixtures)) {
    values <- fixtures[[name]]
    formatted <- sprintf("%.17e", values)
    formatted[is.na(values)] <- "NaN"
    formatted[is.infinite(values) & values > 0] <- "Inf"
    formatted[is.infinite(values) & values < 0] <- "-Inf"
    writeLines(c(name, as.character(length(values)), formatted), connection)
}
close(connection)

# Provenance is regenerated alongside the numerical fixture, without timestamps.
provenance <- c(
    "Oracle: genuine installed fTrading 3042.79, with all 45 source bodies/formals",
    "verified identical to the bundled upstream R functions before generation.",
    paste("R:", R.version.string),
    paste("timeSeries:", packageVersion("timeSeries")),
    paste("timeDate:", packageVersion("timeDate")),
    paste("fBasics:", packageVersion("fBasics")),
    paste("Records:", length(fixtures)),
    paste("Numeric values:", sum(lengths(fixtures))),
    "NA and NaN are both encoded as NaN; infinities retain their sign.",
    "Matrices are flattened in column-major order, matching Fortran storage.",
    "No random-number generator, dependency clone, source repair, or stub is used.",
    "Known original vector-path errors; genuine timeSeries paths are compared:",
    paste(names(vector_failures), vector_failures, sep = ": ")
)
writeLines(provenance, "test/reference/provenance.txt")
cat(paste(provenance, collapse = "\n"), "\n")
