# Run from package root with binom 1.1-2 installed. No RNG is required.
suppressPackageStartupMessages(library(binom))
stopifnot(as.character(packageVersion('binom')) == '1.1.2')
options(digits = 17)
out <- 'test/data'
dir.create(out, showWarnings = FALSE, recursive = TRUE)
write <- function(x, name) write.table(x, file.path(out, name), row.names = FALSE, col.names = FALSE,
                                     quote = FALSE, sep = ',')
methods <- c('agresti-coull','asymptotic','bayes','cloglog','exact','logit','probit','profile','lrt','prop.test','wilson')
rows <- list()
for(n in c(1, 10, 35, 1000)) {
  for(x in unique(c(0, 1, n%/%3, n%/%2, n-1, n))) for(level in c(.8,.95,.99)) {
    for(method in methods) {
      ci <- suppressWarnings(binom.confint(x,n,level,method))
      rows[[length(rows)+1]] <- data.frame(method,x,n,level,ci$mean,ci$lower,ci$upper)
    }
  }
}
write(do.call(rbind,rows),'intervals.csv')
rows <- list()
for(method in methods) for(p in c(0,.01,.37,.5,.91,1)) {
  cov <- binom.coverage(p,10,.95,method)$coverage
  len <- binom.length(p,10,.95,method)$length
  rows[[length(rows)+1]] <- data.frame(method,p,10,.95,cov,len)
}
write(do.call(rbind,rows),'summary.csv')
rows <- list()
for(method in c('cloglog','logit','probit','asymp','lrt','exact')) for(alt in c('two.sided','greater','less')) {
  for(pa in c(.15,.4,.65,.9)) {
    value <- binom.power(pa,30,.4,.05,1.3,alt,method)
    rows[[length(rows)+1]] <- data.frame(method,alt,pa,30,.4,.05,1.3,value)
  }
}
write(do.call(rbind,rows),'power.csv')
rows <- list()
for(type in c('highest','central')) for(x in c(0,2,5,10)) {
  ci <- binom.bayes(x,10,.9,type,1.3,2.1)
  rows[[length(rows)+1]] <- data.frame(type,x,10,.9,1.3,2.1,ci$mean,ci$lower,ci$upper,ci$sig)
}
write(do.call(rbind,rows),'bayes.csv')
rows <- list()
for(x in c(0,1,3,5,10)) for(adj in c(FALSE,TRUE)) for(b in c(FALSE,TRUE)) {
  ci <- binom.lrt(x,10,.95,bayes=b,conf.adj=adj)
  rows[[length(rows)+1]] <- data.frame(x,10,.95,as.integer(b),as.integer(adj),ci$lower,ci$upper)
}
write(do.call(rbind,rows),'lrt.csv')
rows <- list()
for(x in c(0,1,3,5,10)) for(b in c(FALSE,TRUE)) {
  ci <- binom.profile(x,10,.9,maxsteps=100,del=.15,bayes=b)
  rows[[length(rows)+1]] <- data.frame(x,10,.9,as.integer(b),ci$lower,ci$upper)
}
write(do.call(rbind,rows),'profile.csv')
rows <- list()
for(alt in c('two.sided','greater','less')) {
  pa <- if(alt=='less') .3 else .7
  for(exact in c(FALSE,TRUE)) {
    d <- cloglog.sample.size(pa,p=.4,power=.8,alpha=.05,alternative=alt,exact.n=exact,recompute.power=TRUE,phi=1.3)
    rows[[length(rows)+1]] <- data.frame(alt,as.integer(exact),pa,d$n,d$power)
  }
}
write(do.call(rbind,rows),'design.csv')
rows <- list()
for(alt in c('two.sided','greater','less')) {
 d <- cloglog.sample.size(n=100,p=.4,power=.4,alpha=.05,alternative=alt,phi=1.3)
 rows[[length(rows)+1]] <- data.frame(alt,d$p.alt)
}
write(do.call(rbind,rows),'detectable.csv')
u <- seq(.025,.975,by=.05)
counts <- qbinom(u,10,.4)
write(data.frame(u,counts),'simulation_counts.csv')
rows <- list()
for(method in methods) {
 ci <- binom.confint(counts,10,.95,method)
 captured <- sum(.4>=ci$lower & .4<=ci$upper)
 capture <- binom.confint(captured,length(counts),.95,method)
 rows[[length(rows)+1]] <- data.frame(method,captured,capture$mean,capture$lower,capture$upper,mean(ci$upper-ci$lower))
}
write(do.call(rbind,rows),'simulation.csv')
beta_grid <- expand.grid(p=c(1e-15,1e-8,.025,.5,.975,1-1e-8),
                         a=c(.1,.5,1,2,100,10000), b=c(.1,.5,2,100,10000))
beta_grid$q <- qbeta(beta_grid$p,beta_grid$a,beta_grid$b)
write(beta_grid,'beta.csv')
cat('Generated pinned-R reference rows: ', sum(sapply(list.files(out,full.names=TRUE), function(f)length(readLines(f)))), '\n')
