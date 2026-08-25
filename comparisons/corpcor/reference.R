args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else "r_results.csv"
library(corpcor)
source(file.path("..", "common", "read_asset_prices_binary.R"))
prices <- read_asset_prices_binary(file.path("..", "..", "asset_class_etf_prices.bin"))
rets <- apply(prices$prices,2L,function(x) diff(log(x)))
wcheck <- function(x) sum(as.numeric(x)*seq_along(x))
rows <- list()
add <- function(name,fun,reps=30L,atol=1e-9,rtol=1e-8) {
  tm <- system.time(for(i in seq_len(reps))value <- as.numeric(fun()))[["elapsed"]]
  rows[[length(rows)+1L]] <<- data.frame(case=name,value=value,seconds=tm,
    abs_tol=atol,rel_tol=rtol)
}
add("etf_cor_shrink_fixed",function()wcheck(cor.shrink(rets,lambda=.2,verbose=FALSE)))
add("etf_cov_shrink_fixed",function()wcheck(cov.shrink(rets,lambda=.2,
  lambda.var=.1,verbose=FALSE)))
add("etf_invcov_shrink_fixed",function()wcheck(invcov.shrink(rets,lambda=.2,
  lambda.var=.1,verbose=FALSE)))
add("etf_pcor_shrink_fixed",function()wcheck(pcor.shrink(rets,lambda=.2,verbose=FALSE)))
add("etf_cov_shrink_estimated",function()wcheck(cov.shrink(rets,verbose=FALSE)),reps=10L)
add("etf_cov_lambda_estimated",function()attr(cov.shrink(rets,verbose=FALSE),"lambda"),reps=10L)
add("etf_var_lambda_estimated",function()attr(cov.shrink(rets,verbose=FALSE),"lambda.var"),reps=10L)
write.csv(do.call(rbind,rows),out,row.names=FALSE,quote=FALSE)
