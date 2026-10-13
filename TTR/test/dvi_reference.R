# SPDX-License-Identifier: GPL-2.0-or-later
# Regression for sliding-sum rounding changing upstream asymmetric rank ties.
# Run from the package root: Rscript test/dvi_reference.R
suppressPackageStartupMessages(library(TTR))
stopifnot(as.character(packageVersion('TTR')) == '0.24.4')
i <- 1:512
x <- 80 + .033*i + 2.8*sin(.173*i) + cos(.039*i) + .25*cos(.91*i)
z <- DVI(x,n=21,wts=c(.3,.7),smooth=4,magnitude=c(7,19,4),
         stretch=c(5,13,3),exact.multiplier=.3)
con <- file('test/dvi_reference.dat','w')
for (j in seq_along(x)) {
 row <- c(x[j],z[j,])
 writeLines(paste(ifelse(is.na(row),'NaN',sprintf('%.17g',row)),collapse=' '),con)
}
close(con)
