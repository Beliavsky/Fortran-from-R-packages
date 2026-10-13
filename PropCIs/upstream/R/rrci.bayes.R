rrci.bayes <- function(x1,n1,x2,n2,a,b,c,d,conf.level=0.95, nsim = 10000000)
{
  risk.app<- function(a1,b1,c1,d1,conf.level,nsim=nsim)
  {
    z1 <- rbeta(nsim, a1,b1)
    z2 <- rbeta(nsim, c1,d1)
    z <- z1/z2
    z <- sort(z)
    lq <- nsim * (1-conf.level)/2
    uq <- nsim * (1 - (1-conf.level)/2)
    ci <- array(0,2)
    ci[1] <- z[lq]
    ci[2] <- z[uq]
    return(ci)
  }
  
  # Bayes tail interval with beta priors
  
  fct.F1<- function(x,t,a1,b1,a2,b2){
    dbeta(x,a2,b2)*pbeta(x*t,a1,b1)}
  
  fct.F2<- function(x,t,a1,b1,a2,b2){
    dbeta(x,a2,b2)*(1-pbeta(x*t,a1,b1))}
  
  risk.F <- function(t,a1,b1,a2,b2)
  {
    if((0<t) && (t<=1)){
      return(integrate(fct.F1,0,1,t=t,a1=a1,b1=b1,a2=a2,b2=b2)$value)
    }
    else{
      return(1-integrate(fct.F2,0,1/t,t=t,a1=a1,b1=b1,a2=a2,b2=b2)$value)
    }
  }
  
  risk.fct <- function(ab,a1,b1,c1,d1,conf.level)
  {
    abs(risk.F(ab[2],a1,b1,c1,d1) - (1 - (1-conf.level)/2)) +
      abs(risk.F(ab[1],a1,b1,c1,d1) -(1-conf.level)/2)
  }
  

    a1 <- a + x1
    b1 <- b + n1 - x1
    c1 <- c + x2
    d1 <- d + n2 - x2
    start <- risk.app(a1,b1,c1,d1,conf.level, nsim)
    tailci <- optim(start,risk.fct,a1=a1,b1=b1,c1=c1,d1=d1,
                    conf.level=conf.level,control=list(maxit=20000))$par
    if(tailci[1] < 0) tailci[1]  <- 0
    return(tailci)
  }