#' EP distribution
#' @export
#' @name EPd
#' @param x,q vector of quantiles.
#' @param lambda,beta are parameters.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#'  to be the number required.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P\left[ X\leq x\right]}, otherwise, \eqn{P\left[ X>x\right] }.
#' @description
#' Density, distribution function, quantile function and random generation for
#' the EP distribution.
#' @return \code{dEPd} gives the density, \code{pEPd} gives the distribution
#' function, \code{qEPd} gives the quantile function and \code{rEPd} generates
#' random deviates.
#' @details
#' The EP distribution with parameters \eqn{\lambda} and \eqn{\beta},
#' has density
#'  \deqn{f\left( x\right) =\frac{\lambda \beta }
#' {\left( 1-e^{-\lambda }\right) } e^{-\lambda -\beta x+\lambda e^{-\beta x}},}
#' where
#'  \deqn{x>\mathbb{R}_{+},~\beta ,\lambda \in \mathbb{R}_{+}.}
#' @references  Kuş, C., 2007,
#' *A new lifetime distribution*, Computational Statistics & Data Analysis,
#' 51 (9), 4497-4509.
#' @examples
#' library(new.dist)
#' dEPd(1, lambda=2, beta=3)
dEPd<-function(x,lambda,beta,log=FALSE)
{
  if(any(lambda<=0)) {stop("lambda must be > 0")}
  if(any(beta<=0)) {stop("beta must be > 0")}
  enuzun<-max(length(x),length(lambda),length(beta))
  x<-rep(x,enuzun/length(x)+1)[1:enuzun]
  lambda<-rep(lambda, enuzun/length(lambda)+1)[1:enuzun]
  beta<-rep(beta,enuzun/length(beta)+1)[1:enuzun]
  pdf<-NULL
  for (i in 1:enuzun)
  {
    if(x[i]<=0) {pdf[i]<-0} else
    {pdf[i]<-(lambda[i]*beta[i])/(1-exp(-lambda[i]))*exp((-lambda[i])-beta[i]*
                                          x[i]+lambda[i]*exp(-beta[i]*x[i]))}
  }
  if(log==TRUE) pdf<-log(pdf)
  return(pdf)
}
#' EP distribution
#' @export
#' @rdname EPd
#' @examples
#' pEPd(1,lambda=2,beta=3)
pEPd<-function(q,lambda,beta,lower.tail=TRUE,log.p=FALSE)
{
  if(any(lambda<=0)) {stop("lambda must be > 0")}
  if(any(beta<=0)) {stop("beta must be > 0")}
  enuzun<-max(length(q),length(lambda),length(beta))
  q<-rep(q,enuzun/length(q)+1)[1:enuzun]
  lambda<-rep(lambda, enuzun/length(lambda)+1)[1:enuzun]
  beta<-rep(beta,enuzun/length(beta)+1)[1:enuzun]
  cdf<-NULL
  for (i in 1:enuzun)
  {
    if(q[i]>0) cdf[i]<-((exp(lambda[i]*exp(-beta[i]*q[i]))-exp(lambda[i]))*
                          (1-exp(lambda[i]))^(-1)) else cdf[i]<-0
  }
  if(lower.tail==FALSE) cdf<-1-cdf
  if(log.p==TRUE) cdf<-log(cdf)
  return(cdf)
}
#' EP distribution
#' @export
#' @rdname EPd
#' @examples
#' qEPd(.8,lambda=2,beta=3)
qEPd<-function(p,lambda,beta,lower.tail=TRUE)
{
  if(any(p<0)|any(p>1)) {stop("p must be between >= 0 and <= 1")}
  if(any(lambda<=0)) {stop("lambda must be > 0")}
  if(any(beta<=0)) {stop("beta must be > 0")}
  enuzun<-max(length(p),length(lambda),length(beta))
  p<-rep(p,enuzun/length(p)+1)[1:enuzun]
  lambda<-rep(lambda,enuzun/length(lambda)+1)[1:enuzun]
  beta<-rep(beta,enuzun/length(beta)+1)[1:enuzun]
  qfonk<-NULL
  for (i in 1:enuzun)
  {
    qfonk[i]<--(log(log((exp(lambda[i])+p[i]-p[i]*exp(lambda[i])))/lambda[i])/
                  beta[i])
  }
  if(lower.tail==FALSE)
  {
    qfonk[i]<--(log(log((exp(lambda[i])+(1-p[i])-(1-p[i])*exp(lambda[i])))/
                      lambda[i])/beta[i])
  }
  return(qfonk)
}
#' EP distribution
#' @export
#' @rdname EPd
#' @examples
#' rEPd(10,lambda=2,beta=3)
rEPd<-function(n,lambda,beta)
{
  n<-floor(n)
  if(any(n<1)) {stop("n must be >= 1")}
  if(any(lambda<=0)) {stop("lambda must be > 0")}
  if(any(beta<=0)) {stop("beta must be > 0")}
  rn<-qEPd(stats::runif(n),lambda,beta)
  return(rn)
}
