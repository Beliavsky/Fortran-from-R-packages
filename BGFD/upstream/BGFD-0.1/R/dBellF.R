
dBellF<-function (x, a, b, lambda, log = FALSE)
{
G=(1-((1+(x/a)^b))^(-1))
g=((1+(x/a)^b))^(-2) *b*x^(b-1)*a^(-b)


  pdf <- x
  pdf[log == FALSE] <-lambda*g*exp(lambda*(1-G))*exp(-exp(lambda)*(1-exp(-lambda*G)))/(1-(exp(-exp(lambda)+1)))
  pdf[log == TRUE] <-log(lambda)+log(g)+(lambda*(1-G))-exp(lambda)*(1-exp(-lambda*G))-log((1-(exp(-exp(lambda)+1))))
  return(pdf)
}

