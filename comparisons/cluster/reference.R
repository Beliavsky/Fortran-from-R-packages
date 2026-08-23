args <- commandArgs(trailingOnly=TRUE); out <- if(length(args)) args[1] else "r_results.csv"
library(cluster)
x <- cbind(c(0,.2,-.1,.1,5,5.2,4.9,5.1),c(0,-.1,.2,.1,5,4.8,5.2,5.1))
labels <- c(1,1,1,1,2,2,2,2); rows <- list()
add <- function(name,fun,reps=300,atol=1e-9,rtol=1e-8) {
 tm <- system.time(for(j in 1:reps) value <- fun())[["elapsed"]]
 rows[[length(rows)+1]] <<- data.frame(case=name,value=value,seconds=tm,abs_tol=atol,rel_tol=rtol)
}
add("daisy_sum",function() sum(as.matrix(daisy(x))),50000)
add("pam_objective",function() pam(x,2)$objective[2],50000)
add("silhouette_average",function() mean(silhouette(labels,dist(x))[,3]),50000)
add("agnes_coefficient",function() agnes(x,method="average")$ac,50000,1e-8)
add("diana_coefficient",function() diana(x)$dc,50000,1e-8)
write.csv(do.call(rbind,rows),out,row.names=FALSE,quote=FALSE)
