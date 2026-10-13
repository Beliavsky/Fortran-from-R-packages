# Requires the unmodified CRAN PropCIs 0.3-0 package. Run from package root.
options(digits=17)
library(PropCIs)
stopifnot(as.character(packageVersion('PropCIs')) == '0.3.0')
rows <- list()
add <- function(name, args, tol=1e-7) {
  val <- tryCatch(do.call(get(name, asNamespace('PropCIs')), args), error=function(e) e)
  if (inherits(val,'error')) {
    rows[[length(rows)+1]] <<- c(name,paste(args,collapse=','),'ERROR','ERROR',tol,conditionMessage(val))
    return()
  }
  if(is.list(val)) val <- val$conf.int
  if(length(val)==1) val <- rep(val,2)
  rows[[length(rows)+1]] <<- c(name,paste(args,collapse=','),sprintf('%.17g',val[1]),sprintf('%.17g',val[2]),tol,'')
}
for (nm in c('add4ci','addz2ci','exactci','blakerci','midPci','scoreci')) {
  for(x in c(0,1,5,19,20)) add(nm,list(x,20,0.95))
  add(nm,list(1,1,0.9))
  add(nm,list(1,100,0.99))
}
for (nm in c('diffpropci.mp','diffpropci.Wald.mp','scoreci.mp')) {
  for(v in list(c(0,0,20),c(3,8,20),c(0,20,20),c(20,0,20),c(1,0,1),c(0,1,1))) add(nm,as.list(c(v,0.95)))
}
for(v in list(c(3,8),c(0,8),c(8,0),c(1,1))) add('oddsratioci.mp',as.list(c(v,0.95)))
for(nm in c('wald2ci','diffscoreci','orscoreci','riskscoreci')) {
  for(v in list(c(10,40,18,50),c(0,10,0,13),c(10,10,13,13),c(0,10,5,13),c(10,10,5,13),
               c(5,10,0,13),c(5,10,13,13),c(0,10,13,13),c(10,100,90,100),c(90,100,10,100))) {
    if(nm=='wald2ci') { add(nm,c(as.list(c(v,0.95)),list('AC'))); add(nm,c(as.list(c(v,0.95)),list('Wald'))) } else add(nm,as.list(c(v,0.95)),if(nm=='riskscoreci') 2e-4 else 2e-7)
  }
}
for(p in c(0,0.1,0.5,0.9,1)) for(x in c(0,2,10)) add('acceptbin',list(x,10,p))
for(v in list(c(0.25,40,0.36,50,-0.2),c(0,10,0.4,10,0),c(1,10,0,10,0.9),c(0.5,20,0.5,20,0))) add('z2stat',as.list(v),1e-6)
add('limit',list(10,40,18,50,0.95,0.4,0),1e-8)
add('limit',list(10,40,18,50,0.95,0.7,1),1e-8)
set.seed(4281)
for(nm in c('diffci.bayes','rrci.bayes','orci.bayes')) add(nm,list(5,20,10,30,1,1,1,1,0.95,10000),0.002)
add('diffci.bayes.hpd',list(5,20,10,30,1,1,1,1,0.95),0.003)
add('orci.bayes',list(5,20,30,30,1,1,1,1,0.95,10000))
write.table(do.call(rbind,rows),'test/r_reference.tsv',sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
# Independent R distribution definitions, no simulation or package optimization.
beta_cdf <- function(t,s,kind) {
  loc <- digamma(s[3])-digamma(s[4])
  scale <- sqrt(trigamma(s[3])+trigamma(s[4]))
  integrate(function(u) {
    z <- loc+scale*u
    x <- plogis(z)
    weight <- exp(s[3]*plogis(z,log.p=TRUE)+s[4]*plogis(z,lower.tail=FALSE,log.p=TRUE)-lbeta(s[3],s[4])+log(scale))
    target <- switch(kind,x+t,x*t,plogis(z+log(t)))
    weight*pbeta(target,s[1],s[2])
  },-Inf,Inf,rel.tol=1e-10,subdivisions=1000)$value
}
ind <- list()
for(s in list(c(1,1,1,1),c(0.5,0.5,0.5,0.5),c(2.5,11.5,7.2,4.1),c(6,16,31,1),c(20,2,1,9),c(0.2,4,3,0.3))) {
  for(kind in 1:3) {
    interval <- if(kind==1) c(-1,1) else c(-100,100)
    q <- function(p) uniroot(function(q) beta_cdf(if(kind==1) q else exp(q),s,kind)-p,interval,tol=1e-10)$root
    v <- c(q(0.025),q(0.975))
    if(kind!=1) v <- exp(v)
    ind[[length(ind)+1]] <- c(kind,s,v)
  }
}
write.table(do.call(rbind,ind),'test/independent_reference.tsv',sep='\t',row.names=FALSE,col.names=FALSE,quote=FALSE)
