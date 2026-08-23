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
    reps <- if(op == "quantile" && nm %in% c("norm","snorm")) 2000 else if(op == "quantile") 50 else 1000
    tm <- system.time(for(i in seq_len(reps)) value <- wcheck(fun()))[["elapsed"]]
    rows[[length(rows)+1]] <- data.frame(case=paste(nm,op,sep="_"),value=value,
      seconds=tm,abs_tol=5e-5,rel_tol=2e-7)
  }
}
write.csv(do.call(rbind,rows),out,row.names=FALSE,quote=FALSE)
