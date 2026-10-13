# Regenerate committed scalar-reference records with matrixTests 0.2.3.1.
stopifnot(as.character(packageVersion('matrixTests')) == '0.2.3.1')
library(matrixTests)
X <- rbind(c(1.1,3.2,2.5,5.8,4.7,7.4,6.6,8.3,10.2,9.5,12.4,11.6),
           c(0,0,1,1,2,2,3,3,4,5,6,6),
           c(NA,2,3,4,1,NA,2,5,6,7,8,NA),
           1e7+c(1,4,2,8,3,9,5,7,12,6,10,11),
           rep(2,12))
Y <- rbind(c(1.5,2.8,3.5,4.1,5.9,6.3,8.2,7.7,9.8,12.2,10.3,11.1),
           c(1,1,2,3,4,4,5,5,6,7,8,9),
           c(1,2,NA,4,2,3,4,NA,7,9,10,11),
           1e7+c(2,5,1,7,4,10,6,8,11,7,9,13),
           rep(3,12))
g <- rep(1:3, each=4)
t <- seq(0,22,2)
functions <- sub('export\\((.*)\\)', '\\1', readLines('upstream/NAMESPACE'))
functions <- functions[grepl('^row_', functions)]
con <- file('test/reference.txt','w')
fields <- c('obs','obs.x','obs.y','obs.tot','obs.paired','obs.groups','mean','mean.x','mean.y','mean.diff',
'var','var.x','var.y','var.diff','var.pooled','var.ratio','stderr','df','df.num','df.denom','df.between','df.within',
'df.model','df.residual','statistic','pvalue','conf.low','conf.high','cor','skewness','kurtosis',
'sumsq.between','sumsq.within','meansq.between','meansq.within','mesor','amplitude','acrophase','rsquared','period','exact','corrected')
for (f in functions) {
  formal <- names(formals(get(f)))
  modes <- if ('alternative' %in% formal) 1:3 else 1
  if ('exact' %in% formal) modes <- c(modes,4)
  for (mode in modes) {
    args <- list(x=X)
    if ('y' %in% formal) args$y <- Y
    if ('g' %in% formal) args$g <- g
    if ('t' %in% formal) args$t <- t
    if (mode==2) args$alternative <- 'less'
    if (mode==3) args$alternative <- 'greater'
    if (mode %in% 2:3 && 'conf.level' %in% formal) args$conf.level <- .9
    if (mode %in% 2:3 && 'null' %in% formal) args$null <- if(grepl('f_var',f)) 1.3 else .4
    if (mode==4) args$exact <- FALSE
    if (mode==4 && 'correct' %in% formal) args$correct <- FALSE
    result <- suppressWarnings(do.call(f,args))
    for (i in seq_len(nrow(result))) for (field in intersect(fields,names(result))) {
      v <- as.numeric(result[i,field])
      writeLines(paste(f,mode,i,field,format(v,digits=17,scientific=TRUE)),con)
    }
  }
}
close(con)
writeLines(paste('R',getRversion(),'matrixTests',packageVersion('matrixTests'),
                 'matrixStats',packageVersion('matrixStats')),'test/REFERENCE_VERSION.txt')
