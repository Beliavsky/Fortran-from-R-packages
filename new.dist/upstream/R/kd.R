#' Kumaraswamy Distribution
#' @export
#' @name kd
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#'  to be the number required.
#' @param alpha,lambda are non-negative shape parameters.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P\left[ X\leq x\right]}, otherwise, \eqn{P\left[ X>x\right] }.
#' @description
#' Density, distribution function, quantile function and random generation for
#' Kumaraswamy distribution with \code{shape} parameters.
#' @return \code{dkd} gives the density, \code{pkd} gives the distribution
#' function, \code{qkd} gives the quantile function and \code{rkd} generates
#'  random deviates.
#' @details
#' Kumaraswamy distribution with non-negative shape
#' parameters \eqn{\alpha} and \eqn{\lambda} has density
#' \deqn{f\left( x\right) =\alpha \lambda x^{\lambda -1}\left( 1-x^{\lambda }
#' \right)^{\alpha -1},}
#' where
#' \deqn{0<x<1,~~\alpha ,\lambda >0.}
#' @references  Kohansal, A. and Bakouch, H. S., 2021,
#' *Estimation procedures for Kumaraswamy distribution parameters under
#' adaptive type-II hybrid progressive censoring*, Communications in
#' Statistics-Simulation and Computation, 50 (12), 4059-4078.
#' @examples
#' library("new.dist")
#' dkd(0.1,lambda=2,alpha=3)
dkd<-function(x,lambda,alpha,log=FALSE)
{
  if(any(lambda<=0)) {stop("lambda must be > 0")}
  if(any(alpha<=0)) {stop("alpha must be > 0")}
  enuzun <- max(length(x),length(lambda),length(alpha))
  x<-rep(x,enuzun/length(x)+1)[1:enuzun]
  lambda<-rep(lambda, enuzun/length(lambda)+1)[1:enuzun]
  alpha<-rep(alpha,enuzun/length(alpha)+1)[1:enuzun]
  pdf<-NULL
  for (i in 1:enuzun)
  {
    if(x[i]<=0 || x[i]>=1) {pdf[i]<-0} else
    pdf[i]<-alpha[i]*lambda[i]*(x[i]^(lambda[i]-1))*
        (1-x[i]^lambda[i])^(alpha[i]-1)
  }
  if(log==TRUE) pdf<-log(pdf)
  return(pdf)
}
#' Kumaraswamy Distribution
#' @export
#' @rdname kd
#' @examples
#' pkd(0.5,lambda=2,alpha=3)
pkd<-function(q,lambda,alpha,lower.tail=TRUE,log.p=FALSE)
  {
  if(any(lambda<=0)) {stop("lambda must be > 0")}
  if(any(alpha<=0)) {stop("alpha must be > 0")}
    enuzun <- max(length(q),length(lambda),length(alpha))
    q<-rep(q,enuzun/length(q)+1)[1:enuzun]
    lambda<-rep(lambda, enuzun/length(lambda)+1)[1:enuzun]
    alpha<-rep(alpha,enuzun/length(alpha)+1)[1:enuzun]
    cdf<-NULL
    for (i in 1:enuzun)suppressWarnings(
    {
      if (q[i]>0 && q[i]<1) cdf[i]<-1-(1-q[i]^lambda[i])^alpha[i] else cdf[i]<-0
    })
    if(lower.tail==FALSE) cdf<-1-cdf
    if(log.p==TRUE) cdf<-log(cdf)
    return(cdf)
  }
#' Kumaraswamy Distribution
#' @export
#' @rdname kd
#' @examples
#' qkd(.8,lambda=2,alpha=3)
qkd<-function(p,lambda,alpha,lower.tail=TRUE) # 0<p<1, lambda,alpha>0
{
  if(any(p<0)|any(p>1)) {stop("p must be between >= 0 and <= 1")}
  if(any(lambda<=0)) {stop("lambda must be > 0")}
  if(any(alpha<=0)) {stop("alpha must be > 0")}
  enuzun <- max(length(p),length(lambda),length(alpha))
  p<-rep(p,enuzun/length(p)+1)[1:enuzun]
  lambda<-rep(lambda, enuzun/length(lambda)+1)[1:enuzun]
  alpha<-rep(alpha,enuzun/length(alpha)+1)[1:enuzun]
  qfonk<-NULL
  for (i in 1:enuzun)
  {
    qfonk[i]<-exp(log(-exp(log(1-p[i])/alpha[i])+1)/lambda[i])
  }
  if(lower.tail==FALSE)
  {
    qfonk[i]<-exp(log(-exp(log(1-(1-p[i]))/alpha[i])+1)/lambda[i])
  }
  return(qfonk)
}
#' Kumaraswamy Distribution
#' @export
#' @rdname kd
#' @examples
#' rkd(10,lambda=2,alpha=3)
rkd<-function(n,lambda,alpha)
  {
    n<-floor(n)
    if(any(n<1)) {stop("n must be >= 1")}
    if(any(lambda<=0)) {stop("lambda must be > 0")}
    if(any(alpha<=0)) {stop("alpha must be > 0")}
    suppressWarnings({
    rn<-qkd(stats::runif(n),lambda,alpha)})
    return(rn)
  }
