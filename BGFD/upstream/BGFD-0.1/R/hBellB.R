
hBellB<-function (x, a, b, k, lambda, log = FALSE,log.p = FALSE, lower.tail = TRUE) 
{
G=(1-((1+(x/a)^b))^(-k))
g=k*((1+(x/a)^b))^(-k-1) *b*x^(b-1)*a^(-b)

  pdf <- x
  pdf[log == FALSE] <-lambda*g*exp(lambda*(1-G))*exp(-exp(lambda)*(1-exp(-lambda*G)))/(1-(exp(-exp(lambda)+1)))
  pdf[log == TRUE] <-log(lambda)+log(g)+(lambda*(1-G))-exp(lambda)*(1-exp(-lambda*G))-log((1-(exp(-exp(lambda)+1))))
  
  cdf <- x
  cdf[log.p == FALSE & lower.tail == TRUE] <- (1-exp(-exp(lambda)*(1-exp(-lambda*G))))/(1-(exp(-exp(lambda)+1))) 
  cdf[log.p == TRUE & lower.tail == TRUE] <- log(1-exp(-exp(lambda)*(1-exp(-lambda*G))))-log(1-(exp(-exp(lambda)+1)))  
  cdf[log.p == FALSE & lower.tail == FALSE] <- ((exp(-exp(lambda)*(1-exp(-lambda*G))))-(exp(-exp(lambda)+1)))/(1-(exp(-exp(lambda)+1)))
  cdf[log.p == TRUE & lower.tail == FALSE] <- log((exp(-exp(lambda)*(1-exp(-lambda*G))))-(exp(-exp(lambda)+1)))-log(1-(exp(-exp(lambda)+1)))
  
  hrf<-pdf/(1-cdf)
  return(hrf)
}