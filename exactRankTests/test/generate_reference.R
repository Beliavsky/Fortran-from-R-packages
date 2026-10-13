# Reference generator, exactRankTests 0.8-37. Run from the package directory.
# SPDX-License-Identifier: GPL-2.0-or-later
suppressPackageStartupMessages(library(exactRankTests))
stopifnot(packageDescription('exactRankTests')$Version == '0.8-37')
con <- file('test/reference_values.dat', 'w')
emit <- function(x) {
  txt <- format(as.numeric(x), digits=17, scientific=TRUE, trim=TRUE)
  txt[is.na(x)] <- 'NaN'
  writeLines(paste(txt, collapse=' '), con)
}
x <- c(1,2,2,3,5)
y <- c(0,2,4,4)
a <- c(-2,-1,0,1,3,3)
for (ex in c(TRUE,FALSE)) for (alt in c('two.sided','less','greater')) {
  z <- wilcox.exact(x,y,alternative=alt,exact=ex)
  emit(c(z$statistic,z$p.value,if(ex) z$pointprob else 0))
  z <- wilcox.exact(a,alternative=alt,exact=ex)
  emit(c(z$statistic,z$p.value,if(ex) z$pointprob else 0))
  z <- ansari.exact(x,y,alternative=alt,exact=ex)
  emit(c(z$statistic,z$p.value))
  z <- perm.test(x,y,alternative=alt,exact=ex)
  emit(c(z$statistic,z$p.value,if(ex) z$pointprob else 0))
  z <- perm.test(a,alternative=alt,exact=ex)
  emit(c(z$statistic,z$p.value,if(ex) z$pointprob else 0))
}
for(alt in c('two.sided','less','greater')) {
  z <- wilcox.exact(x,y,alternative=alt,exact=TRUE,conf.int=TRUE,conf.level=.8)
  emit(c(z$conf.int,z$estimate))
  z <- wilcox.exact(a,alternative=alt,exact=TRUE,conf.int=TRUE,conf.level=.8)
  emit(c(z$conf.int,z$estimate))
}
for (ex in c(TRUE,FALSE)) for (tol in list(NULL,0.01)) {
  z <- perm.test(c(-1.2,.5,2.7,4), c(.1,.3,2.3), exact=ex,tol=tol)
  emit(c(z$statistic,z$p.value))
  z <- perm.test(c(-2.1,-1.2,0,.7,3.3),exact=ex,tol=tol)
  emit(c(z$statistic,z$p.value))
}
y <- c(2,1,2,4,6)
for (int in c(FALSE,TRUE)) for (type in c('Data','Wilcoxon','NormalQuantile','AnsariBradley','Median','Savage','ConSal')) emit(cscores(y,type=type,int=int))
emit(cscores.Surv(cbind(c(1,2,2,3,5),c(1,0,1,1,0))))
emit(cscores.Surv(cbind(c(1,2,2,3,5),c(1,0,1,1,0)),int=TRUE))
scores <- c(-1,0,.5,2.5,2.5)
support <- sort(unique(combn(scores,2,sum)))
emit(support)
emit(dperm(support,scores,2))
for (alt in c('less','greater','two.sided')) emit(pperm(support,scores,2,alternative=alt))
emit(qperm(c(0,.1,.25,.5,.75,.99,1),scores,2))
for (alt in c('two.sided','less','greater')) {
  z <- suppressWarnings(ansari.exact(c(1,2,3,4,5),c(1,2,4,6,9),alternative=alt,exact=TRUE,conf.int=TRUE,conf.level=.8))
  emit(c(z$conf.int,z$estimate))
}
z <- wilcox.exact(c(-1.5,-.5,-2.5,2,-.5,4,3.5,0,-3,-.5,-3.5),mu=1,exact=TRUE,conf.int=TRUE,conf.level=.9)
emit(c(z$conf.int,z$estimate))
close(con)
