diffci.bayes.hpd <- function(x1,n1,x2,n2,a,b,c,d,conf.level=0.95)
{
  diff.app.hpd <- function(a1,b1,c1,d1,conf.level)
  {
    mu1 <- a1/(a1+b1)
    mu2 <- c1/(c1+d1)
    v1 <- mu1*(1-mu1)/(a1 + b1)
    v2 <- mu2*(1-mu2)/(c1 + d1)
    z  <- qnorm(1-(1-conf.level)/2)
    ci <- array(0,2)
    ci[1] <- mu1-mu2 - z*sqrt(v1+v2)
    ci[2] <- mu1-mu2 + z*sqrt(v1+v2)
    return(ci)
  }
  
  # Bayes HPD interval with beta priors
  fct.F1<- function(x,t,a1,b1,c1,d1){
    dbeta(x,c1,d1)*pbeta(x+t,a1,b1)}
  
  fct.F2<- function(x,t,a1,b1,c1,d1){
    dbeta(x,c1,d1)*(1-pbeta(x+t,a1,b1))}
  
  diff.F <- function(t,a1,b1,c1,d1)
  { if(t < 0)
    Fvalue <- integrate(fct.F1,-t,1,t=t,a1=a1,b1=b1,c1=c1,d1=d1)$value
  else
    Fvalue <- 1-integrate(fct.F2,0,1-t,t=t,a1=a1,b1=b1,c1=c1,d1=d1)$value
  return(Fvalue)
  }
  
  fct.f<- function(x,t,a1,b1,c1,d1){
    dbeta(x,c1,d1)*dbeta(x+t,a1,b1)
  }
  
  diff.f <- function(t,a1,b1,c1,d1)
  {
    if(t < -1) fvalue <- 100
    else if(t > 1)  fvalue <- 100
    else if((t >= -1) && (t <= 0))
      fvalue <- integrate(fct.f,-t,1,t=t,a1=a1,b1=b1,c1=c1,d1=d1)$value
    else
      fvalue <- integrate(fct.f,0,1-t,t=t,a1=a1,b1=b1,c1=c1,d1=d1)$value
    return(fvalue)
  }
  
  
  diff <- function(ab,a1,b1,c1,d1,conf.level)
  {
    1000*abs(diff.F(ab[2],a1,b1,c1,d1) -
               diff.F(ab[1],a1,b1,c1,d1) - conf.level)+
      abs(diff.f(ab[1],a1,b1,c1,d1) -
            diff.f(ab[2],a1,b1,c1,d1))
  }
  
  y <- x1
  if( y > n1/2 ) {
    x2 <- n2-x2
    x1 <- n1-x1
  }
  a1 <- a + x1
  b1 <- b + n1 - x1
  c1 <- c + x2
  d1 <- d + n2 - x2
  start <- diff.app.hpd(a1,b1,c1,d1,conf.level)
  hdrci <- optim(start,diff,a1=a1,b1=b1,c1=c1,d1=d1,conf.level=conf.level)$par
  if(hdrci[1] < -1) hdrci[1]  <- -1
  if(hdrci[2] >  1) hdrci[1]  <- 1
  if( y > n1/2 ) {
    ci <- hdrci
    hdrci[1] <- -ci[2]
    hdrci[2] <- -ci[1]
  }
  return(hdrci)
}