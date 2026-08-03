#'@title Data-driven selection of p  for the test of randomness
#'
#'@description This function uses a AIC/BIC type criterion to select p based on the data.
#'
#'@param X           Time series
#'@param p0          Minimum value of p (default is 2)
#'@param d           Maximum value of p (default is 5)
#'@param q           Constant for selecting between AIC and BIC type penalty (default is 2.4)
#'@param lambda      Penalty term (default is 0.25); small values lead to p=d, large value lead to p=p0
#'
#'
#'@return \item{p}{Selected value of p}
#'
#'@references B.R Nasri (2021). Tests of serial dependence for arbitrary distributions
#'

#'@examples
#' X <- SimAR1Poisson(c(5,0.2),100)
#' out <- select_p(X)
#'
#'@export
#'
#'


select_p = function(X,p0=2,d=5,q=2.4,lambda=0.25){
  n = length(X)
  fn=preparedata(X)$fn
  m = (1-sum(fn^2))/6

  out= TestIndSerCopula(X,d,d,1)
  card = out$card
  cvm = out$stat$cvm /(m^card)
  k = length(card)
  M = max(cvm)

  if(M <= q*log(n))
    {cte = lambda*log(n)}
  else
    { cte =2*lambda}


  S=sum(cvm)

  L=d+1-p0
  V= numeric(L)

  V[L] = S-k*cte

  for( i in 1:(L-1))
    {
       p  = p0+(i-1);
      out = TestIndSerCopula(X,p,p,1)
     card = out$card
     cvm  = out$stat$cvm /(m^card)
      k   = length(card)
     S    = sum(cvm)
    V[i]  = S- k*cte
  }


  ind = sort(V,index.return=TRUE)

  p = ind$ix[L]+ p0-1
  return(p)

}





