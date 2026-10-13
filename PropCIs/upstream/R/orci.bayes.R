orci.bayes <- function(x1,n1,x2,n2,a,b,c,d,conf.level=0.95, nsim = 10000000)
{
  or.app<- function(a1,b1,c1,d1,conf.level,nsim=nsim)
  {
    z1 <- rf(nsim, 2*a1,2*b1)
    z2 <- rf(nsim, 2*c1,2*d1)
    a <- (d1/c1)/(b1/a1)
    z <- a*z1/z2
    z <- sort(z)
    lq <- nsim * (1-conf.level)/2
    uq <- nsim * (1 - (1-conf.level)/2)
    ci <- array(0,2)
    ci[1] <- z[lq]
    ci[2] <- z[uq]
    return(ci)
  }
  
  
  # Bayes tail interval with beta priors
  
  fct.F<- function(x,t,a1,b1,a2,b2){
    c <- (b2/a2)/(b1/a1)
    df(x,2*a2,2*b2)*pf(x*t/c,2*a1,2*b1)
  }
  
  or.F <- function(t,a1,b1,a2,b2)
  {
    return(integrate(fct.F,0,Inf,t=t,a1=a1,b1=b1,a2=a2,b2=b2)$value)
  }
  
  or.fct <- function(ab,a1,b1,c1,d1,conf.level)
  {
    abs(or.F(ab[2],a1,b1,c1,d1) - (1 - (1-conf.level)/2))+
      abs(or.F(ab[1],a1,b1,c1,d1) - (1-conf.level)/2)
  }
  

    if(x2!=n2){
      a1 <- a + x1
      b1 <- b + n1 - x1
      c1 <- c + x2
      d1 <- d + n2 - x2
      start <- or.app(a1,b1,c1,d1,conf.level,nsim)
      tailci <- optim(start,or.fct,a1=a1,b1=b1,c1=c1,d1=d1,
                      conf.level=conf.level,control=list(maxit=20000))$par
      if(tailci[1] < 0) tailci[1]  <- 0 }
    else{
      a1 <- a + n1 - x1
      b1 <- b +  x1
      c1 <- c + n2 - x2
      d1 <- d + x2
      start <- or.app(a1,b1,c1,d1,conf.level,nsim)
      tailci1 <- optim(start,or.fct,a1=a1,b1=b1,c1=c1,d1=d1,
                       conf.level=conf.level,control=list(maxit=20000))$par
      if(tailci[1] < 0) tailci[1]  <- 0
      tailci <- array(0,2)
      tailci[1] <- 1/ tailci1[2]
      tailci[2] <- 1/ tailci1[1]
    }
    return(tailci)
  }