args <- commandArgs(trailingOnly=TRUE); out <- if(length(args)) args[1] else "r_results.csv"
library(fracdiff)
n <- 4096; i <- 1:n; x <- sin(.017*i)+.35*cos(.0031*i)+.0002*i
wcheck <- function(z) sum(z*seq_along(z)); rows <- list()
add <- function(name,fun,reps=10,atol=1e-7,rtol=1e-8) {
 tm <- system.time(for(j in 1:reps) value <- fun())[["elapsed"]]
 rows[[length(rows)+1]] <<- data.frame(case=name,value=value,seconds=tm,abs_tol=atol,rel_tol=rtol)
}
add("diffseries_d025",function() wcheck(diffseries(x,.25)),10,1e-5)
add("diffseries_d070",function() wcheck(diffseries(x,.70)),10,1e-5)
add("gph_d",function() fdGPH(x)$d,30)
add("sperio_d",function() fdSperio(x)$d,30)
write.csv(do.call(rbind,rows),out,row.names=FALSE,quote=FALSE)
