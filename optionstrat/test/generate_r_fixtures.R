# SPDX-License-Identifier: GPL-3.0-only
# Reproduce fixtures with genuine CRAN optionstrat 1.4.1.
# Run from the package root: Rscript test/generate_r_fixtures.R <output-directory>
library(optionstrat)
stopifnot(as.character(packageVersion('optionstrat')) == '1.4.1')
out <- commandArgs(trailingOnly = TRUE)
if (length(out) != 1) stop('Provide one output directory')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_numeric <- function(x, name) {
  lines <- apply(x, 1, function(row) paste(sprintf('%.17g', row), collapse = ','))
  writeLines(lines, file.path(out, name))
}
inputs <- rbind(c(100,100,.2,.25,.03,.01), c(100,95,.3,2,-.01,.02),
                c(90,110,.4,.6,.02,.005), c(120,90,.15,1,.05,.03),
                c(25,30,.7,.03,.02,0), c(10,10,1.2,4,.03,.01),
                c(100,110,.2,0,.02,0), c(100,90,.2,0,.02,0))
analytics <- t(apply(inputs, 1, function(a) {
  args <- as.list(a)
  c(a, unlist(do.call(calleval,args)), unlist(do.call(puteval,args)))
}))
write_numeric(analytics, 'analytics.csv')
extra <- t(apply(inputs[1:6,], 1, function(a) {
  args <- as.list(a)
  c(do.call(lambda,c(list('call'),args)), do.call(lambda,c(list('put'),args)),
    iv.calc('call',do.call(callpremium,args),a[1],a[2],a[4],a[5],a[6]),
    iv.calc('put',do.call(putpremium,args),a[1],a[2],a[4],a[5],a[6]))
}))
write_numeric(extra, 'extra.csv')
probinputs <- rbind(c(100,90,110,.2,0,30,.04,262,.2),
                   c(100,80,115,0,.01,45,0,262,.75),
                   c(50,30,70,.3,.02,20,-.02,250,.9),
                   c(100,90,110,.2,0,0,0,262,.75))
probabilities <- t(apply(probinputs, 1, function(a) {
  args <- list(spot=a[1],asd=a[4],dsd=a[5],dte=a[6],mean=a[7],tradedays=a[8])
  quantiles <- c(args,list(p=a[9],quantile=TRUE))
  c(a,do.call(prob.below,c(args,list(upper=a[3]))),
    do.call(prob.above,c(args,list(lower=a[2]))),
    do.call(prob.btwn,c(args,list(lower=a[2],upper=a[3]))),
    do.call(prob.below,quantiles)$price,do.call(prob.above,quantiles)$price,
    do.call(prob.btwn,quantiles)$price)
}))
write_numeric(probabilities, 'probabilities.csv')
spreads <- rbind(unlist(vertical('call',100,90,110,.25,.02,.2,sigma2=.25,vol=.3,d=.01)),
                 unlist(vertical('call',100,110,90,.25,.02,.2,sigma2=.25,vol=.3,d=.01)),
                 unlist(vertical('put',100,90,110,.25,.02,.2,sigma2=.25,vol=.3,d=.01)),
                 unlist(vertical('put',100,110,90,.25,.02,.2,sigma2=.25,vol=.3,d=.01)),
                 unlist(dv(100,80,90,110,120,.25,.02,.2,.21,.22,.23,vol=.3,d=.01)))
write_numeric(spreads, 'spreads.csv')
writeLines(c(paste('R',getRversion()),paste('optionstrat',packageVersion('optionstrat'))),file.path(out,'versions.txt'))
