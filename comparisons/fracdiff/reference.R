args <- commandArgs(trailingOnly=TRUE); out <- if(length(args)) args[1] else "r_results.csv"
library(fracdiff)
n <- 4096; i <- 1:n; x <- sin(.017*i)+.35*cos(.0031*i)+.0002*i
wcheck <- function(z) sum(z*seq_along(z)); rows <- list()
add <- function(name,fun,reps=10,atol=1e-7,rtol=1e-8) {
 tm <- system.time(for(j in 1:reps) value <- fun())[["elapsed"]]
 rows[[length(rows)+1]] <<- data.frame(case=name,value=value,seconds=tm,abs_tol=atol,rel_tol=rtol)
}
add("diffseries_d025",function() wcheck(diffseries(x,.25)),100,1e-5)
add("diffseries_d070",function() wcheck(diffseries(x,.70)),100,1e-5)
add("gph_d",function() fdGPH(x)$d,30)
add("sperio_d",function() fdSperio(x)$d,30)

source(file.path("..", "common", "read_asset_prices_binary.R"))
prices <- read_asset_prices_binary(file.path("..", "..", "asset_class_etf_prices.bin"))
spy <- diff(log(prices$prices[,"SPY"]))
spy_squared <- (spy - mean(spy))^2
add("etf_spy_diffseries_d025", function() wcheck(diffseries(spy,.25)), 30, 1e-7)
add("etf_spy_gph_d", function() fdGPH(spy)$d, 10)
add("etf_spy_sperio_d", function() fdSperio(spy)$d, 10)
add("etf_spy_squared_gph_d", function() fdGPH(spy_squared)$d, 10)
add("etf_spy_squared_sperio_d", function() fdSperio(spy_squared)$d, 10)
write.csv(do.call(rbind,rows),out,row.names=FALSE,quote=FALSE)
