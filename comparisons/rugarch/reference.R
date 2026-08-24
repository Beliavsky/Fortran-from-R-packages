args <- commandArgs(trailingOnly=TRUE)
out <- if(length(args)) args[1] else "r_results.csv"
library(rugarch)
rows <- list(); x <- seq(-4,4,length.out=1001); p <- seq(.001,.999,length.out=1001)
wcheck <- function(z) sum(z*seq_along(z))
specs <- list(norm=c(5,1), std=c(7,1), ged=c(1.6,1), snorm=c(5,1.3),
              sstd=c(7,1.3), sged=c(1.6,1.3), jsu=c(1.8,.5))
for(nm in names(specs)) {
  shape <- specs[[nm]][1]; skew <- specs[[nm]][2]
  for(op in c("density","cdf","quantile")) {
    fun <- switch(op,
      density=function() ddist(nm,x,shape=shape,skew=skew),
      cdf=function() pdist(nm,x,shape=shape,skew=skew),
      quantile=function() qdist(nm,p,shape=shape,skew=skew))
    reps <- if(op == "quantile" && nm %in% c("norm","snorm")) 2000 else if(op == "quantile") 500 else 1000
    tm <- system.time(for(i in seq_len(reps)) value <- wcheck(fun()))[["elapsed"]]
    rows[[length(rows)+1]] <- data.frame(case=paste(nm,op,sep="_"),value=value,
      seconds=tm,abs_tol=5e-5,rel_tol=2e-7)
  }
}

source(file.path("..", "common", "read_asset_prices_binary.R"))
prices <- read_asset_prices_binary(file.path("..", "..", "asset_class_etf_prices.bin"))
spy <- 100*diff(log(prices$prices[,"SPY"]))
fixed_filter <- function(model, gamma1=NULL) {
  pars <- list(mu=.04, omega=.02, alpha1=.08, beta1=.90)
  if (!is.null(gamma1)) pars$gamma1 <- gamma1
  spec <- ugarchspec(
    variance.model=list(model=model, garchOrder=c(1,1)),
    mean.model=list(armaOrder=c(0,0), include.mean=TRUE),
    distribution.model="norm",
    fixed.pars=pars)
  ugarchfilter(spec, data=spy, filter.control=list(rec.init="all"))
}
add_filter_case <- function(name, model, gamma1=NULL, value_fun, reps=20,
                            atol=1e-6, rtol=1e-8) {
  tm <- system.time(for (j in seq_len(reps)) {
    filtered <- fixed_filter(model, gamma1)
    value <- value_fun(filtered)
  })[["elapsed"]]
  rows[[length(rows)+1]] <<- data.frame(case=name,value=value,seconds=tm,
    abs_tol=atol,rel_tol=rtol)
}
add_filter_case("etf_spy_sgarch_sigma", "sGARCH", value_fun=function(x) wcheck(as.numeric(sigma(x))))
conditional_loglik <- function(x) {
  vol <- as.numeric(sigma(x))
  eps <- as.numeric(residuals(x))
  sum(dnorm(eps[-1L]/vol[-1L], log=TRUE)-log(vol[-1L]))
}
add_filter_case("etf_spy_sgarch_loglik", "sGARCH", value_fun=conditional_loglik,
                rtol=5e-8)
add_filter_case("etf_spy_gjrgarch_sigma", "gjrGARCH", gamma1=.04,
                value_fun=function(x) wcheck(as.numeric(sigma(x))))
add_filter_case("etf_spy_gjrgarch_loglik", "gjrGARCH", gamma1=.04,
                value_fun=conditional_loglik, rtol=5e-8)
write.csv(do.call(rbind,rows),out,row.names=FALSE,quote=FALSE)
