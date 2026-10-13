# Recreate pinned posterior 1.7.0 oracle data: Rscript test/reference_diagnostics.R
# Run from package root. R and posterior are test-generation dependencies only.
stopifnot(as.character(packageVersion("posterior")) == "1.7.0")
suppressPackageStartupMessages(library(posterior))
set.seed(481516)
base <- matrix(rnorm(320), 80, 4)
ar <- base
anti <- base
for (j in 1:4) for (i in 2:80) {
  ar[i,j] <- 0.95 * ar[i-1,j] + base[i,j]
  anti[i,j] <- -0.85 * anti[i-1,j] + base[i,j]
}
cases <- list(independent=base, positive_ar=ar, negative_ar=anti,
              shifted=sweep(base,2,c(0,0,2,5),`+`), scaled=sweep(base,2,c(1,1,3,10),`*`),
              ties=round(base), odd=base[1:79,], onechain=base[,1,drop=FALSE],
              constant=matrix(3,10,4), different_constants=matrix(rep(1:4,each=10),10),
              one_constant=cbind(rep(2,80),base[,2:4]), single_iteration=matrix(1:4,1),
              short2=base[1:2,], squared_constant=matrix(rep(c(-1,1),40),20),
              short3=base[1:3,],short4=base[1:4,],short5=base[1:5,],short6=base[1:6,],
              missing=replace(base,13,NA), infinite=replace(base,13,Inf))
con <- file("test/reference_diagnostics.dat", "w")
writeLines(as.character(length(cases)),con)
number <- function(x) ifelse(is.na(x),"NaN",ifelse(is.infinite(x),ifelse(x>0,"Inf","-Inf"),sprintf("%.17g",x)))
for (nm in names(cases)) {
  x <- cases[[nm]]
  expect <- suppressWarnings(c(rhat_basic(x),rhat_basic(x,split=FALSE),rhat(x),
    ess_basic(x),ess_basic(x,split=FALSE),ess_bulk(x),ess_tail(x),ess_mean(x),ess_median(x),ess_sd(x),
    mcse_mean(x),mcse_median(x),mcse_sd(x),ess_quantile(x,c(0,.3,1)),mcse_quantile(x,c(0,.3,1)),
    rhat_nested(x,rep(seq_len(ncol(x)),length.out=ncol(x)))))
  writeLines(nm,con)
  writeLines(paste(nrow(x),ncol(x)),con)
  writeLines(paste(number(c(x)),collapse=" "),con)
  writeLines(paste(number(expect),collapse=" "),con)
}
close(con)
