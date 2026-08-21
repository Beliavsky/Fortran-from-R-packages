
dBellW<-function (x, alpha,beta,lambda, log = FALSE)
{
  G=(1-exp(-alpha*x^beta))
g=alpha*beta*x^(beta-1)*exp(-alpha*x^beta)
  pdf <- x
  pdf[log == FALSE] <-lambda*g*exp(lambda*(1-G))*exp(-exp(lambda)*(1-exp(-lambda*G)))/(1-(exp(-exp(lambda)+1)))
  pdf[log == TRUE] <-log(lambda)+log(g)+(lambda*(1-G))-exp(lambda)*(1-exp(-lambda*G))-log((1-(exp(-exp(lambda)+1))))
  return(pdf)
}

