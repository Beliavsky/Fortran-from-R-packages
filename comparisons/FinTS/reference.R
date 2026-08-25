args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else "r_results.csv"
suppressPackageStartupMessages(library(FinTS))
source(file.path("..", "common", "read_asset_prices_binary.R"))
prices <- read_asset_prices_binary(file.path("..", "..", "asset_class_etf_prices.bin"))
spy <- diff(log(prices$prices[,"SPY"]))
spy2 <- (spy-mean(spy))^2
wcheck <- function(x) sum(as.numeric(x)*seq_along(x))
rows <- list()
add <- function(name, fun, reps=100L, atol=1e-9, rtol=1e-8) {
  tm <- system.time(for(i in seq_len(reps)) value <- as.numeric(fun()))[["elapsed"]]
  rows[[length(rows)+1L]] <<- data.frame(case=name,value=value,seconds=tm,
    abs_tol=atol,rel_tol=rtol)
}
add("etf_spy_acf20", function() wcheck(Acf(spy,lag.max=20,plot=FALSE)$acf))
add("etf_spy_squared_acf20", function() wcheck(Acf(spy2,lag.max=20,plot=FALSE)$acf))
add("etf_spy_pacf20", function() wcheck(Acf(spy,lag.max=20,type="partial",plot=FALSE)$acf))
add("etf_spy_ljung_box_stat", function() unname(AutocorTest(spy,lag=20)$statistic))
add("etf_spy_ljung_box_pvalue", function() AutocorTest(spy,lag=20)$p.value)
add("etf_spy_arch12_stat", function() unname(ArchTest(spy,lags=12)$statistic), reps=20L)
add("etf_spy_arch12_pvalue", function() ArchTest(spy,lags=12)$p.value, reps=20L)
write.csv(do.call(rbind,rows),out,row.names=FALSE,quote=FALSE)
