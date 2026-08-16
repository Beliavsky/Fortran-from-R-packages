#' Discrete Lindley Distribution
#' @export
#' @name dLd2
#' @param x,q vector of quantiles.
#' @param theta a parameter.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#'  to be the number required.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P\left[ X\leq x\right]}, otherwise, \eqn{P\left[ X>x\right] }.
#' @description
#' Density, distribution function, quantile function and random generation for
#' the discrete Lindley distribution.
#' @return \code{ddLd2} gives the density, \code{pdLd2} gives the distribution
#' function, \code{qdLd2} gives the quantile function and \code{rdLd2} generates
#' random deviates.
#' @details
#' the discrete Lindley distribution with a parameter \eqn{\theta},
#' has density
#'  \deqn{f\left( x\right) =\frac{\lambda ^{x}}{1+\theta }
#'  \left( \theta \left(1-2\lambda \right) +\left( 1-\lambda \right)
#'  \left( 1+\theta x\right)\right),}
#' where
#'  \deqn{x=0,1,2,...~,\lambda =\exp \left( -\theta \right) ~and~\theta >0.}
#' @references  Bakouch, H. S., Jazi, M. A. and Nadarajah, S., 2014,
#' *A new discrete distribution, Statistics*, 48 (1), 200-240.
#' @examples
#' library(new.dist)
#' ddLd2(2,theta=2)
ddLd2<-function(x,theta,log=FALSE)
{
  x<-floor(x)
  if(any(theta<=0)) {stop("theta must be > 0")}
  enuzun <- max(length(x),length(theta))
  x<-rep(x,enuzun/length(x)+1)[1:enuzun]
  theta<-rep(theta,enuzun/length(theta)+1)[1:enuzun]
  pdf<-NULL
  lambda<-exp(-theta)
  for (i in 1:enuzun)
  {
    if(x[i]<0) {pdf[i]<-0} else
    {pdf[i]<-(lambda[i]^x[i]/(1+theta[i]))*(theta[i]*(1-2*lambda[i])+
                                              (1-lambda[i])*(1+theta[i]*x[i]))}
  }
  if(log==TRUE) pdf<-log(pdf)
  return(pdf)
}
#' Discrete Lindley Distribution
#' @export
#' @rdname dLd2
#' @examples
#' pdLd2(1,theta=2)
pdLd2<-function(q,theta,lower.tail=TRUE,log.p=FALSE)
{
  q<-floor(q)+1
  if(any(theta<=0)) {stop("theta must be > 0")}
  enuzun <- max(length(q),length(theta))
  q<-rep(q,enuzun/length(q)+1)[1:enuzun]
  theta<-rep(theta,enuzun/length(theta)+1)[1:enuzun]
  cdf<-NULL
  lambda<-exp(-theta)
  for (i in 1:enuzun)
  {

    if(q[i]>=0) cdf[i]<-1-((1+theta[i]+theta[i]*(q[i]))/(1+theta[i]))*
        lambda[i]^(q[i]) else cdf[i]<-0
  }
  if(lower.tail==FALSE) cdf<-1-cdf
  if(log.p==TRUE) cdf<-log(cdf)
  return(cdf)
}
#' Discrete Lindley Distribution
#' @export
#' @rdname dLd2
#' @examples
#' qdLd2(.5,theta=2)
qdLd2<-function(p,theta,lower.tail=TRUE)
{
  if(any(p<0)|any(p>1)) {stop("p must be between >= 0 and <= 1")}
  if(any(theta<=0)) {stop("theta must be > 0")}
  enuzun <- max(length(p),length(theta))
  p<-rep(p,enuzun/length(p)+1)[1:enuzun]
  theta<-rep(theta,enuzun/length(theta)+1)[1:enuzun]
  quant<-NULL
  for (i in 1:enuzun)
  {
    x<-0
    t<-0
    while(t<p[i]){

      t<-pdLd2(x,theta[i])
      x<-x+1
    }
    quant[i]<-x-1
  }
  if(lower.tail==FALSE)
  {
    for (i in 1:enuzun)
    {
      x<-0
      t<-0
      while(t<(1-p[i])){
        t<-pdLd2(x,theta[i])
        x<-x+1
      }
      quant[i]<-x-1
    }
  }
  return(quant)
}
#' Discrete Lindley Distribution
#' @export
#' @rdname dLd2
#' @examples
#' rdLd2(10,theta=1)
rdLd2<-function(n,theta)
  {
    n<-floor(n)
    if(any(n<1)) {stop("n must be >= 1")}
    if(any(theta<=0)) {stop("theta must be > 0")}
    suppressWarnings({
    rn<-qdLd2(stats::runif(n),theta)})
    return(rn)
  }
