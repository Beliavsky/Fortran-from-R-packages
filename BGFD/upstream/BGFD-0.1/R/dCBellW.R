
dCBellW<-function (x, alpha,beta,lambda, log = FALSE)
{
  G=(1-exp(-alpha*x^beta))
g=alpha*beta*x^(beta-1)*exp(-alpha*x^beta)
  pdf <- x
  pdf[log == FALSE] <- lambda*g*exp(lambda*G)*exp(exp(lambda*G)-1)/(exp(exp(lambda)-1)-1)
  pdf[log == TRUE] <- log(lambda)+log(g)+(lambda*G)+(exp(lambda*G)-1)-log((exp(exp(lambda)-1)-1))
  return(pdf)
}

