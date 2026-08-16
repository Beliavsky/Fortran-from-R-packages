#' Unit Inverse Gaussian Distribution
#' @export
#' @name uigd
#' @param x,q vector of quantiles.
#' @param mu a mean parameter.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#'  to be the number required.
#' @param lambda a scale parameter.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P\left[ X\leq x\right]}, otherwise, \eqn{P\left[ X>x\right] }.
#' @description
#' Density, distribution function, quantile function and random generation for
#' the Unit Inverse Gaussian distribution \code{mean} and \code{scale}.
#' @return \code{duigd} gives the density, \code{puigd} gives the distribution
#' function, \code{quigd} gives the quantile function and \code{ruigd} generates
#'  random deviates.
#' @details
#' The Unit Inverse Gaussian distribution \code{scale}
#' parameter \eqn{\lambda} and \code{mean}
#' parameter \eqn{\mu}, has density
#' \deqn{f\left( x\right) =\sqrt{\frac{\lambda }{2\pi }}
#' \frac{1}{x^{3/2}}e^{-\frac{ \lambda }{2\mu ^{2}x}\left( x-\mu \right) ^{2}},}
#' where
#' \deqn{x>0,~\mu ,\lambda >0.}
#' @references  Ghitany, M., Mazucheli, J., Menezes, A. and Alqallaf, F., 2019,
#' *The unit-inverse Gaussian distribution: A new alternative to two-parameter*
#' *distributions on the unit interval, Communications in Statistics-Theory and*
#'  *Methods*, 48 (14), 3423-3438.
#' @examples
#' library(new.dist)
#' duigd(1, mu=2, lambda=3)
duigd<-function(x,mu,lambda=1,log=FALSE)
{
  if(any(mu<=0)) {stop("mu must be > 0")}
  if(any(lambda<=0)) {stop("lambda must be > 0")}
  enuzun <- max(length(x),length(mu),length(lambda))
  x<-rep(x,enuzun/length(x)+1)[1:enuzun]
  mu<-rep(mu, enuzun/length(mu)+1)[1:enuzun]
  lambda<-rep(lambda,enuzun/length(lambda)+1)[1:enuzun]
  pdf<-NULL
  for (i in 1:enuzun)
  {
    if(x[i]<=0) {pdf[i]<-0} else
    pdf[i]<-((lambda[i]/(2*pi))^(1/2))*(1/(x[i]^(3/2)))*exp(-(lambda[i]/
                                          (2*mu[i]^2*x[i]))*(x[i]-mu[i])^2)
  }
  if(log==TRUE) pdf<-log(pdf)
  return(pdf)
}
#' Unit Inverse Gaussian Distribution
#' @export
#' @rdname uigd
#' @examples
#' puigd(1,mu=2,lambda=3)
puigd<-function(q,mu,lambda=1,lower.tail=TRUE,log.p=FALSE)
{
  if(any(mu<=0)) {stop("mu must be > 0")}
  if(any(lambda<=0)) {stop("lambda must be > 0")}
  enuzun <- max(length(q),length(mu),length(lambda))
  q<-rep(q,enuzun/length(q)+1)[1:enuzun]
  mu<-rep(mu, enuzun/length(mu)+1)[1:enuzun]
  lambda<-rep(lambda,enuzun/length(lambda)+1)[1:enuzun]
  cdf<-NULL
  for (i in 1:enuzun)
  {
    if (q[i]>0) (cdf[i]<-stats::pnorm((lambda[i]/q[i])^(1/2)*(q[i]/mu[i]-1))
    +exp(2*lambda[i]/mu[i])*stats::pnorm(-(lambda[i]/q[i])^(1/2)*
                                           (q[i]/mu[i]+1))) else cdf[i]<-0
  }
  if(lower.tail==FALSE) cdf<-1-cdf
  if(log.p==TRUE) cdf<-log(cdf)
  return(cdf)
}
#' Unit Inverse Gaussian Distribution
#' @export
#' @rdname uigd
#' @examples
#' quigd(.1,mu=2,lambda=3)
quigd<-function(p,mu,lambda=1,lower.tail=TRUE)
{
  if(any(p<0)|any(p>1)) {stop("p must be between >= 0 and <= 1")}
  if(any(mu<=0)) {stop("mu must be > 0")}
  if(any(lambda<=0)) {stop("lambda must be > 0")}
  enuzun <- max(length(p),length(mu),length(lambda))
  p<-rep(p,enuzun/length(p)+1)[1:enuzun]
  mu<-rep(mu, enuzun/length(mu)+1)[1:enuzun]
  lambda<-rep(lambda,enuzun/length(lambda)+1)[1:enuzun]
  kok<-NULL
  for (i in 1:enuzun)
  {
    if (lower.tail==TRUE) {
    Y<-function(x)
    {
      stats::pnorm((lambda[i]/x)^(1/2)*(x/mu[i]-1))+exp(2*lambda[i]/mu[i])*
        stats::pnorm(-(lambda[i]/x)^(1/2)*(x/mu[i]+1))-p[i]
    }}
    else
    {
      Y<-function(x)
      {
        stats::pnorm((lambda[i]/x)^(1/2)*(x/mu[i]-1))+exp(2*lambda[i]/mu[i])*
          stats::pnorm(-(lambda[i]/x)^(1/2)*(x/mu[i]+1))-(1-p[i])
      }
    }
    kok[i]<-(stats::uniroot(Y,c(0,100000)))$root
  }
  return(kok)
}
#' Unit Inverse Gaussian Distribution
#' @export
#' @rdname uigd
#' @examples
#' ruigd(10,mu=2,lambda=3)
ruigd<-function(n,mu,lambda=1)
{
  n<-floor(n)
  if(any(n<1)) {stop("n must be >= 1")}
  if(any(mu<=0)) {stop("mu must be > 0")}
  if(any(lambda<=0)) {stop("lambda must be > 0")}
  rn<-quigd(stats::runif(n),mu,lambda)
  return(rn)
}
