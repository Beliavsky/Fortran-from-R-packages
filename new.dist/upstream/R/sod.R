#' Standard Omega Distribution
#' @export
#' @name sod
#' @param x,q vector of quantiles.
#' @param alpha,beta are parameters.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#'  to be the number required.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P\left[ X\leq x\right]}, otherwise, \eqn{P\left[ X>x\right] }.
#' @description
#' Density, distribution function, quantile function and random generation for
#' the Standard Omega distribution.
#' @return \code{dsod} gives the density, \code{psod} gives the distribution
#' function, \code{qsod} gives the quantile function and \code{rsod} generates
#' random deviates.
#' @details
#' The Standard Omega distribution with parameters
#' \eqn{\alpha} and \eqn{\beta}, has density
#' \deqn{f\left( x\right) =\alpha \beta x^{\beta -1}\frac{1}{1-x^{2\beta }}
#' \left( \frac{1+x^{\beta }}{1-x^{\beta }}\right) ^{-\alpha /2},}
#' where
#' \deqn{0<x<1,~\alpha ,\beta >0.}
#' @references  Birbiçer, İ. and Genç, A. İ., 2022,
#' *On parameter estimation of the standard omega distribution*. Journal of
#' Applied Statistics, 1-17.
#' @examples
#' library(new.dist)
#' dsod(0.4, alpha=1, beta=2)
dsod<-function(x,alpha,beta,log=FALSE)
  {
    if(any(alpha<=0)) {stop("alpha must be > 0")}
    if(any(beta<=0)) {stop("beta must be > 0")}
    enuzun<-max(length(x),length(alpha),length(beta))
    x<-rep(x,enuzun/length(x)+1)[1:enuzun]
    alpha<-rep(alpha, enuzun/length(alpha)+1)[1:enuzun]
    beta<-rep(beta,enuzun/length(beta)+1)[1:enuzun]
    pdf<-NULL
    for (i in 1:enuzun)
    {
      suppressWarnings(
      if(x[i]<=0 | x[i]>=1) {pdf[i]<-0} else
      {pdf[i]<-(alpha[i]*beta[i]*x[i]^(beta[i]-1))*(1/(1-x[i]^(2*beta[i])))*
        ((1+x[i]^beta[i])/(1-x[i]^beta[i]))^(-alpha[i]/2)})
    }
    if(log==TRUE) pdf<-log(pdf)
    return(pdf)
}
#' Standard Omega Distribution
#' @export
#' @rdname sod
#' @examples
#' psod(0.4, alpha=1, beta=2)
psod<-function(q,alpha,beta,lower.tail=TRUE,log.p=FALSE)
{
  if(any(alpha<=0)) {stop("alpha must be > 0")}
  if(any(beta<=0)) {stop("beta must be > 0")}
  enuzun <- max(length(q),length(alpha),length(beta))
  q<-rep(q,enuzun/length(q)+1)[1:enuzun]
  alpha<-rep(alpha, enuzun/length(alpha)+1)[1:enuzun]
  beta<-rep(beta,enuzun/length(beta)+1)[1:enuzun]
  cdf<-NULL
  for (i in 1:enuzun)
  {
    if (q[i] > 0 && q[i] < 1) {cdf[i]<-1-((1+q[i]^beta[i])/
                              (1-q[i]^beta[i]))^(-alpha[i]/2)} else (cdf[i] <-0)
  }
  if(lower.tail==FALSE) cdf<-1-cdf
  if(log.p==TRUE) cdf<-log(cdf)
  return(cdf)
}
#' Standard Omega Distribution
#' @export
#' @rdname sod
#' @examples
#' qsod(.8, alpha=1, beta=2)
qsod<-function(p,alpha,beta,lower.tail=TRUE)
{
  if(any(p<0)|any(p>1)) {stop("p must be between >= 0 and <= 1")}
  if(any(alpha<=0)) {stop("alpha must be > 0")}
  if(any(beta<=0)) {stop("beta must be > 0")}
  enuzun <- max(length(p),length(alpha),length(beta))
  p<-rep(p,enuzun/length(p)+1)[1:enuzun]
  alpha<-rep(alpha, enuzun/length(alpha)+1)[1:enuzun]
  beta<-rep(beta,enuzun/length(beta)+1)[1:enuzun]
  qfonk<-NULL
  for (i in 1:enuzun)
  {
    qfonk[i]<-exp(log((exp(-(2*log(1-p[i])/alpha[i]))-1)/
                        (exp(-(2*log(1-p[i])/alpha[i]))+1))/beta[i])
  }
  if(lower.tail==FALSE)
  {
    qfonk[i]<-exp(log((exp(-(2*log(1-(1-p[i]))/alpha[i]))-1)/
                        (exp(-(2*log(1-(1-p[i]))/alpha[i]))+1))/beta[i])
  }
  return(qfonk)
}
#' Standard Omega Distribution
#' @export
#' @rdname sod
#' @examples
#' rsod(10, alpha=1, beta=2)
rsod<-function(n,alpha,beta)
  {
    n<-floor(n)
    if(any(n<1)) {stop("n must be >= 1")}
    if(any(alpha<=0)) {stop("alpha must be > 0")}
    if(any(beta<=0)) {stop("beta must be > 0")}
    suppressWarnings({
    rn<-qsod(stats::runif(n),alpha,beta)})
    return(rn)
  }
