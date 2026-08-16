#' Ram Awadh Distribution
#' @export
#' @name RA
#' @param x,q vector of quantiles.
#' @param theta a scale parameter.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#'  to be the number required.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P\left[ X\leq x\right]}, otherwise, \eqn{P\left[ X>x\right] }.
#' @description
#' Density, distribution function, quantile function and random generation for
#' a Ram Awadh distribution with parameter \code{scale}.
#' @return \code{dRA} gives the density, \code{pRA} gives the distribution
#'  function, \code{qRA} gives the quantile function and \code{rRA}
#'  generates random deviates.
#' @details
#' Ram Awadh distribution with \code{scale} parameter
#' \eqn{\theta}, has density
#'  \deqn{f\left( x\right) =\frac{\theta ^{6}}{\theta ^{6}+120}
#'  \left( \theta+x^{5}\right) e^{-\theta x},}
#' where
#'  \deqn{x>0,~\theta >0.}
#' @references  Shukla, K. K., Shanker, R. and Tiwari, M. K., 2022,
#' *A new one parameter discrete distribution and its applications*, Journal
#' of Statistics and Management Systems, 25 (1), 269-283.
#' @examples
#' library(new.dist)
#' dRA(1,theta=2)
dRA<-function(x,theta=1,log=FALSE)
{
  if(any(theta<=0)) {stop("theta must be > 0")}
  enuzun <- max(length(x),length(theta))
  x<-rep(x,enuzun/length(x)+1)[1:enuzun]
  theta<-rep(theta,enuzun/length(theta)+1)[1:enuzun]
  pdf<-NULL
  for (i in 1:enuzun)
  {
    if(x[i]<=0) {pdf[i]<-0} else
    {pdf[i]<-(theta[i]^6/(theta[i]^6+120))*(theta[i]+x[i]^5)*
      exp(-theta[i]*x[i])}
  }
  if(log==TRUE) pdf<-log(pdf)
  return(pdf)
}
#' Ram Awadh Distribution
#' @export
#' @rdname RA
#' @examples
#' pRA(1,theta=2)
pRA<-function(q,theta=1,lower.tail=TRUE,log.p=FALSE)
{
  if(any(theta<=0)) {stop("theta must be > 0")}
  enuzun<-max(length(q),length(theta))
  q<-rep(q,enuzun/length(q)+1)[1:enuzun]
  theta<-rep(theta,enuzun/length(theta)+1)[1:enuzun]
  cdf<-NULL
  for (i in 1:enuzun)
  {
    if(q[i]>0) cdf[i]<-1-(1+((theta[i]*q[i]*(theta[i]^4*q[i]^4+5*theta[i]^3*
          q[i]^3+20*theta[i]^2*q[i]^2+60*theta[i]*q[i]+120)/(theta[i]^6+120))))*
        exp(-theta[i]*q[i]) else cdf[i]<-0
  }
  if(lower.tail==FALSE) cdf<-1-cdf
  if(log.p==TRUE) cdf<-log(cdf)
  return(cdf)
}
#' Ram Awadh Distribution
#' @export
#' @rdname RA
#' @examples
#' qRA(.1,theta=1)
qRA<-function(p,theta=1,lower.tail=TRUE)
{
  if(any(p<0)|any(p>1)) {stop("p must be between >= 0 and <= 1")}
  if(any(theta<=0)) {stop("theta must be > 0")}
  enuzun<-max(length(p),length(theta))
  p<-rep(p,enuzun/length(p)+1)[1:enuzun]
  theta<-rep(theta,enuzun/length(theta)+1)[1:enuzun]
  kok<-NULL
  for (i in 1:enuzun)suppressWarnings(
  {
    Ex<-(theta[i]^6+2*3*4*5*6)/(theta[i]*(theta[i]^6+120))
    Y<-function(x)
    {
      abs((1-(1+((theta[i]*x*(theta[i]^4*x^4+5*theta[i]^3*x^3+20*theta[i]^2*
              x^2+60*theta[i]*x+120)/(theta[i]^6+120))))*exp(-theta[i]*x))-p[i])
    }
    if(lower.tail==FALSE)
    {
      Y<-function(x)
      {
        abs((1-(1+((theta[i]*x*(theta[i]^4*x^4+5*theta[i]^3*x^3+20*theta[i]^2*
          x^2+60*theta[i]*x+120)/(theta[i]^6+120))))*exp(-theta[i]*x))-(1-p[i]))
      }
    }
    kok[i]<-stats::optim(Ex,Y)$par
  })
  return(kok)
}
#' Ram Awadh Distribution
#' @export
#' @rdname RA
#' @examples
#' rRA(10,theta=1)
rRA<-function(n,theta=1)
{
  n<-floor(n)
  if(any(n<1)) {stop("n must be >= 1")}
  if(any(theta<=0)) {stop("theta must be > 0")}
  rn<-qRA(stats::runif(n),theta)
  return(rn)
}
