#' Uniform-Geometric Distribution
#' @export
#' @name ugd
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
#' the Uniform-Geometric distribution.
#' @return \code{dugd} gives the density, \code{pugd} gives the distribution
#' function, \code{qugd} gives the quantile function and \code{rugd} generates
#' random deviates.
#' @details
#' The Uniform-Geometric distribution with shape parameter \eqn{\theta}, has
#' density
#' \deqn{f\left( x\right) =\theta \left( 1-\theta \right) ^{x-1}LerchPhi
#' \left[ \left(1-\theta \right) ,1,x\right],}
#' where
#' \deqn{LerchPhi\left( z,a,v\right) =\sum_{n=0}^{\infty }\frac{z^{n}}
#' {\left(v+n\right) ^{a}}}
#' and
#' \deqn{x=1,2,...~,~~0<\theta <1.}
#' @references Akdoğan, Y., Kuş, C., Asgharzadeh, A., Kınacı, İ., and Sharafi,
#'  F. (2016).
#' *Uniform-geometric distribution*. Journal of Statistical Computation and
#' Simulation, 86(9), 1754-1770.
#' @examples
#' library(new.dist)
#' dugd(1, theta=0.5)
dugd<-function(x,theta,log=FALSE)
{
  x<-floor(x)
  if(any(theta<=0)|any(theta>=1)) {stop("theta must be between  0 > and < 1")}
  enuzun <- max(length(x),length(theta))
  x<-rep(x,enuzun/length(x)+1)[1:enuzun]
  theta<-rep(theta,enuzun/length(theta)+1)[1:enuzun]
  pdf<-NULL
  for (i in 1:enuzun)
  {
    if(x[i]<=0) {pdf[i]<-0} else
    {pdf[i]<-theta[i]*(1-theta[i])^(x[i]-1)*VGAM::lerch((1-theta[i]),1,x[i])}
  }
  if(log==TRUE) pdf<-log(pdf)
  return(pdf)
}
#' Uniform-Geometric Distribution
#' @export
#' @rdname ugd
#' @examples
#' pugd(1,theta=.5)
pugd<-function(q,theta,lower.tail=TRUE,log.p=FALSE)
{
  q<-floor(q)
  if(any(theta<=0)|any(theta>=1)) {stop("theta must be between  0 > and < 1")}
  enuzun <- max(length(q),length(theta))
  q<-rep(q,enuzun/length(q)+1)[1:enuzun]
  theta<-rep(theta,enuzun/length(theta)+1)[1:enuzun]
  cdf<-NULL
  for (i in 1:enuzun)
  {
    if(q[i]>=1) cdf[i]<-1-theta[i]*(1-theta[i])^q[i]*((1/theta[i])-q[i]*
                              VGAM::lerch(1-theta[i],1,q[i]+1)) else cdf[i]<-0
  }
  if(lower.tail==FALSE) cdf<-1-cdf
  if(log.p==TRUE) cdf<-log(cdf)
  return(cdf)
}
#' Uniform-Geometric Distribution
#' @export
#' @rdname ugd
#' @examples
#' qugd(0.6,theta=.1)
qugd<-function(p,theta,lower.tail=TRUE)
{
  if(any(p<0)|any(p>1)) {stop("p must be between >= 0 and <= 1")}
  if(any(theta<=0)|any(theta>=1)) {stop("theta must be between  0 > and < 1")}
  enuzun <- max(length(p),length(theta))
  p<-rep(p,enuzun/length(p)+1)[1:enuzun]
  theta<-rep(theta, enuzun/length(theta)+1)[1:enuzun]
  quant<-NULL
  for (i in 1:enuzun)
  {
    x<-0
    t<-0
    while(t<p[i]){
      t<-pugd(x,theta[i])
      x<-x+1
    }
    quant[i]<-x-1
  }
  if(lower.tail==FALSE){
    for (i in 1:enuzun){
      x<-0
      t<-0
      while(t<(1-p[i])){
        t<-pugd(x,theta[i])
        x<-x+1
      }
      quant[i] <- x-1}
  }
  {
    return(quant)
  }
}
#' Uniform-Geometric Distribution
#' @export
#' @rdname ugd
#' @examples
#' rugd(10,theta=.1)
rugd<-function(n,theta)
  {
    n<-floor(n)
    if(any(n<1)) {stop("n must be >= 1")}
    if(any(theta<=0)|any(theta>=1)) {stop("theta must be between  0 > and < 1")}
    suppressWarnings({
    rn<-qugd(stats::runif(n),theta)})
    return(rn)
  }
