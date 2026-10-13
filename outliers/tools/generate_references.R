# Rebuild deterministic fixtures using actual CRAN outliers 0.15 and R stats.
library(outliers)
stopifnot(as.character(packageVersion('outliers')) == '0.15')
args <- commandArgs(TRUE)
root <- if (length(args)) args[1] else '.'
fmt <- function(x) paste(ifelse(is.na(x), 'NaN', sprintf('%.17g', x)), collapse=' ')
safe <- function(expr) suppressWarnings(tryCatch(expr, error=function(e) NaN))
con <- file(file.path(root,'test/distributions.dat'),'w')
record <- function(code,n,type,k,a,b,z) writeLines(fmt(c(code,n,type,k,a,b,z)),con)
pp <- c(0.005,0.01,0.02,0.025,0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95)
for (ty in c(0,10,11,12,20,21,22)) for (n in 3:30) {
  minimum <- switch(as.character(ty),'0'=3,'10'=3,'11'=4,'12'=5,'20'=4,'21'=5,'22'=6)
  if(n<minimum) next
  for(p in unique(c(-.1,0,pp,1,1.1,seq(.015,.985,length.out=19))))
    record(1,n,ty,0,p,0,safe(qdixon(p,n,ty)))
  for(q in unique(c(-.1,0,1,1.1,seq(.007,.993,length.out=29),safe(qdixon(pp,n,ty)))))
    record(2,n,ty,0,q,0,safe(pdixon(q,n,ty)))
}
for (ty in c(10,11,20)) for (n in unique(c(3:30,50,100,1000))) {
  if(ty==20 && (n<4 || n>30)) next
  ps <- unique(c(0,.0001,.001,.01,.025,.05,.1,.15,.2,.4,.6,.8,.9,.95,.975,.99,.999,.9999,.99999,1-1e-12,1))
  for(p in ps) record(3,n,ty,0,p,0,safe(qgrubbs(p,n,ty)))
  qs <- if(ty==20) seq(-.1,1.1,length.out=23) else seq(0,1.01*sqrt(2*(n-1)),length.out=21)
  for(q in unique(c(qs,safe(qgrubbs(ps[ps<1],n,ty)))))
    record(4,n,ty,0,q,0,safe(pgrubbs(q,n,ty)))
}
for(n in c(2,2.5,3,5,10,30,100)) for(k in c(2,3,5,10,30)) {
  for(p in c(0,1e-8,.01,.1,.5,.9,.95,.99,1-1e-10,1))
    record(5,0,0,k,p,n,safe(qcochran(p,n,k)))
  for(q in c(0,1e-8,.01,.1,.2,.5,.9,.99,1))
    record(6,0,0,k,q,n,safe(pcochran(q,n,k)))
}
close(con)

samples <- list(c(1,2,3,4,5,20),c(-40,-20,1,2,3,4,5,6),c(-20,1,2,3,4,5,20),
 c(1,1,1,2,3,20,20),c(NA,1,2,4,7,30,NaN),c(1,2,9),c(1,2,3,10),
 c(1,2),c(3),seq(-2,2,length.out=30)^3,c(-100,seq(-1,1,length.out=14),80,90),rep(2,6))
con <- file(file.path(root,'test/tests.dat'),'w')
for(x in samples) for(opp in c(FALSE,TRUE)) for(two in c(FALSE,TRUE)) {
 for(ty in c(10,11,20)) {
   z <- safe(grubbs.test(x,ty,opp,two))
   if(is.list(z)) {
     g <- if(ty==20) unname(z$statistic['U']) else unname(z$statistic['G'])
     vals <- c(g,unname(z$statistic['U']),z$p.value)
   } else vals <- rep(NaN,3)
   writeLines(fmt(c(1,length(x),ty,opp,two,0,vals)),con)
   writeLines(fmt(x),con)
 }
 for(ty in c(0,10,11,12,20,21,22)) {
   z <- safe(dixon.test(x,ty,opp,two))
   vals <- if(is.list(z)) c(z$statistic,NaN,z$p.value) else rep(NaN,3)
   writeLines(fmt(c(2,length(x),ty,opp,two,0,vals)),con)
   writeLines(fmt(x),con)
 }
 for(v in c(0,3)) {
   z <- if(v==0) safe(chisq.out.test(x,opposite=opp)) else safe(chisq.out.test(x,v,opp))
   vals <- if(is.list(z)) c(z$statistic,NaN,z$p.value) else rep(NaN,3)
   writeLines(fmt(c(3,length(x),0,opp,two,v,vals)),con)
   writeLines(fmt(x),con)
 }
}
close(con)
con <- file(file.path(root,'test/scores.dat'),'w')
for(x in samples) for(ty in seq_along(c('z','t','chisq','iqr','mad'))) {
 for(p in c(NA,0,.1,.5,.93,1)) for(lim in c(NA,0,1.5)) {
   kind <- c('z','t','chisq','iqr','mad')[ty]
   y <- safe(scores(x,kind,prob=p,lim=lim))
   if(length(y)==1 && is.na(y)) y <- rep(NaN,length(x))
   writeLines(fmt(c(length(x),ty,p,lim)),con)
   writeLines(fmt(x),con)
   writeLines(fmt(as.numeric(y)),con)
 }
}
close(con)
cat('Reference fixtures generated with',R.version.string,'and outliers 0.15\n')
