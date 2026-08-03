#'@title Quantile function of margins
#'
#'@description This function computes the quantile of seven cdf used in the simulatuons of Nasri (2022).
#'
#'@param u Vector of probabilities
#'@param k Marginal distribution: [1] Bernoulli(0.8), [2] Poisson(6), [3] Negative binomial with r = 1.5, p = 0.2, [4] Zero-inflated Poisson (10) with w = 0.1 and P(6.67) otherwise, [5] Zero-inflated Gaussian, [6] Discretized Gaussian, [7] Discrete Pareto(1)
#'
#'@return \item{x}{Vector of quantiles}
#'
#'@author Bouchra R. Nasri January 2021
#'@references B.R Nasri (2022). Tests of serial dependence for arbitrary distributions
#'
#'@examples
#' x = Finv(runif(100),2)
#'@export


Finv=  function(u,k)
{

  w = 0.1
  c0 = exp(-6.67)
  p0 = w+(1-w)*c0
  c1 = 0.5*(1-w)



g4 = function(z){return(max(c0,z))}
g5a = function(u){return(min(0.5,u/(1-w)))}
g5b = function(u){return(max(w,(u-w)/(1-w)))}

F1 = function(u){ return(as.double((u>0.2)))}   #Bernoulli(0.8)
F2 = function(u){ return( qpois(u,6)  ) }       # Poisson(6)
F3 = function(u){ return( qnbinom(u,1.5,0.2)) } # Negative binomial with r = 1.5, p = 0.2 (mean = 6)
F4 = function(u){ p = sapply((u-w)/(1-w),g4)
                  return( (u>p0) *qpois(p,6.67)) } # Zero-inflated Poisson (10) with w = 0.1 and P(6.67) otherwise
F5 = function(u){ pa = sapply(u,g5a)
                  pb = sapply(u,g5b)
                 return((u< c1) *qnorm( pa )+(u>w+c1)*qnorm(pb))} # Zero-inflated Gaussian
F6 = function(u){ return( floor(200*qnorm(u) ) )} # Discretized Gaussian
F7 = function(u){ return( -1 +ceiling(1 /(1-u)))}    # Bad luck distribution  = Pareto(1)

#F5 = @(u) (u< c1).*norminv( min(0.5,u/(1-w)))+(u>w+c1).*norminv(max(w,(u-w)/(1-w))) ;% Zero-inflated Gaussian
switch(k,
  { x = F1(u)},
  { x = F2(u)},
  { x = F3(u)},
  { x = F4(u)},
  { x = F5(u)},
  { x = F6(u)},
  { x = F7(u)}
)
return(x)
}
