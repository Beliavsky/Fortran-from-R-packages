#' Two-Parameter Rayleigh Distribution
#' @export
#' @name tprd
#' @param x,q vector of quantiles.
#' @param lambda a scale parameter.
#' @param mu a location parameter.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#'  to be the number required.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P\left[ X\leq x\right]}, otherwise, \eqn{P\left[ X>x\right] }.
#' @description
#' Density, distribution function, quantile function and random generation
#' for the Two-Parameter Rayleigh  distribution with parameters \code{location}
#' and \code{scale}.
#' @return \code{dtprd} gives the density, \code{ptprd} gives the distribution
#' function, \code{qtprd} gives the quantile function and \code{rtprd} generates
#'  random deviates.
#' @details
#' The Two-Parameter Rayleigh distribution with \code{scale} parameter
#' \eqn{\lambda} and \code{location} parameter \eqn{\mu}, has density
#'  \deqn{f\left( x\right) =2\lambda \left( x-\mu \right) e^{-\lambda
#'  \left( x-\mu\right) ^{2}},}
#' where
#'  \deqn{x>\mu ,~\lambda >0.}
#' @references  Dey, S., Dey, T. and Kundu, D., 2014,
#' *Two-parameter Rayleigh distribution: different methods of estimation*,
#' American Journal of Mathematical and Management Sciences, 33 (1), 55-74.
#' @examples
#' library(new.dist)
#' dtprd(5, lambda=4, mu=4)
dtprd<-function(x,lambda=1,mu,log=FALSE)
{
  if(any(lambda<=0)) {stop("lambda must be > 0")}
  enuzun<-max(length(x),length(lambda),length(mu))
  x<-rep(x,enuzun/length(x)+1)[1:enuzun]
  mu<-rep(mu, enuzun/length(mu)+1)[1:enuzun]
  lambda<-rep(lambda,enuzun/length(lambda)+1)[1:enuzun]
  pdf<-NULL
  for (i in 1:enuzun)
  {
    if(x[i]<=mu[i]) {pdf[i]<-0} else
    {pdf[i]<-2*lambda[i]*(x[i]-mu[i])*exp(-lambda[i]*(x[i]-mu[i])^2)}
  }
  if(log==TRUE) pdf<-log(pdf)
  return(pdf)
}
#' Two-Parameter Rayleigh Distribution
#' @export
#' @rdname tprd
#' @examples
#' ptprd(2,lambda=2,mu=1)
ptprd<-function(q,lambda=1,mu,lower.tail=TRUE,log.p=FALSE)
{
  if(any(lambda<=0)) {stop("lambda must be > 0")}
  enuzun<-max(length(q),length(lambda),length(mu))
  q<-rep(q,enuzun/length(q)+1)[1:enuzun]
  lambda<-rep(lambda, enuzun/length(lambda)+1)[1:enuzun]
  mu<-rep(mu,enuzun/length(mu)+1)[1:enuzun]
  cdf<-NULL

  for (i in 1:enuzun)
  {
    if(q[i]>mu[i]) cdf[i]<-1-exp(-lambda[i]*(q[i]-mu[i])^2) else cdf[i]<-0
  }
  if(lower.tail==FALSE) cdf<-1-cdf
  if(log.p==TRUE) cdf<-log(cdf)
  return(cdf)
}
#' Two-Parameter Rayleigh Distribution
#' @export
#' @rdname tprd
#' @examples
#' qtprd(.5,lambda=2,mu=1)
qtprd<-function(p,lambda=1,mu,lower.tail=TRUE)
{
  if(any(p<0)|any(p>1)) {stop("p must be between >= 0 and <= 1")}
  if(any(lambda<=0)) {stop("lambda must be > 0")}
  enuzun<-max(length(p),length(lambda),length(mu))
  p<-rep(p,enuzun/length(p)+1)[1:enuzun]
  lambda<-rep(lambda,enuzun/length(lambda)+1)[1:enuzun]
  mu<-rep(mu,enuzun/length(mu)+1)[1:enuzun]
  qfonk<-NULL
  for(i in 1:enuzun)
  {
    qfonk[i]<- (lambda[i]*mu[i]+(-lambda[i]*log(1-p[i]))^(1/2))/lambda[i]
  }
  if(lower.tail==FALSE)
  {
    qfonk[i]<- (lambda[i]*mu[i]+(-lambda[i]*log(1-(1-p[i])))^(1/2))/lambda[i]
  }
  return(qfonk)
}
#' Two-Parameter Rayleigh Distribution
#' @export
#' @rdname tprd
#' @examples
#' rtprd(10,lambda=2,mu=1)
rtprd<-function(n,lambda=1,mu)
{
  n<-floor(n)
  if(any(n<1)) {stop("n must be >= 1")}
  if(any(lambda<=0)) {stop("lambda must be > 0")}
  rn<-qtprd(stats::runif(n),lambda,mu)
  return(rn)
}
