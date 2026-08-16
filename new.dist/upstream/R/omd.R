#' Muth Distribution
#' @export
#' @name omd
#' @param x,q vector of quantiles.
#' @param alpha a parameter.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#'  to be the number required.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P\left[ X\leq x\right]}, otherwise, \eqn{P\left[ X>x\right] }.
#' @description
#' Density, distribution function, quantile function and random generation for
#' on the Muth distribution.
#' @return \code{domd} gives the density, \code{pomd} gives the distribution
#' function, \code{qomd} gives the quantile function and \code{romd} generates
#' random deviates.
#' @details
#' The Muth distribution with a parameter \eqn{\alpha}, has
#' density
#' \deqn{f\left( x\right) =\left( e^{\alpha x}-
#' \alpha \right) e^{\alpha x-\left(1/\alpha \right) \left( e^{\alpha x}-
#' 1\right) },}
#' where
#' \deqn{x>0,~\alpha \in \left( 0,1\right]. }
#' @references  Jodrá, P., Jiménez-Gamero, M. D. and Alba-Fernández, M. V., 2015,
#' *On the Muth distribution, Mathematical Modelling and Analysis*, 20 (3),
#' 291-310.
#' @examples
#' library(new.dist)
#' domd(1,alpha=.2)
domd<-function(x,alpha,log=FALSE)
{
  if(any(alpha<=0)|any(alpha>1)) {stop("alpha must be > 0 and <= 1")}
  enuzun <- max(length(x),length(alpha))
  x<-rep(x,enuzun/length(x)+1)[1:enuzun]
  alpha<-rep(alpha,enuzun/length(alpha)+1)[1:enuzun]
  pdf<-NULL
  for (i in 1:enuzun)
  {
    if(x[i]<=0) {pdf[i]<-0} else
    {pdf[i]<-(exp(alpha[i]*x[i])-alpha[i])*exp(alpha[i]*x[i]-(1/alpha[i])*
                                                 (exp(alpha[i]*x[i])-1))}
  }
  if(log==TRUE) pdf<-log(pdf)
  return(pdf)
}
#' Muth Distribution
#' @export
#' @rdname omd
#' @examples
#' pomd(1,alpha=.2)
pomd<-function(q,alpha,lower.tail=TRUE,log.p=FALSE)
{
  if(any(alpha<=0)|any(alpha>1)) {stop("alpha must be > 0 and <= 1")}
  enuzun <- max(length(q),length(alpha))
  q<-rep(q,enuzun/length(q)+1)[1:enuzun]
  alpha<-rep(alpha,enuzun/length(alpha)+1)[1:enuzun]
  cdf<-NULL
  for (i in 1:enuzun)
  {
    if(q[i]>0) cdf[i]<-1-exp(alpha[i]*q[i]-(1/alpha[i])*(exp(alpha[i]*q[i])-1))
    else cdf[i]<-0
  }
  if(lower.tail==FALSE) cdf<-1-cdf
  if(log.p==TRUE) cdf<-log(cdf)
  return(cdf)
}
#' Muth Distribution
#' @export
#' @rdname omd
#' @examples
#' qomd(.8,alpha=.1)
qomd<-function(p,alpha,lower.tail=TRUE)
{
  if(any(p<0)|any(p>1)) {stop("p must be between >= 0 and <= 1")}
  if(any(alpha<=0)|any(alpha>1)) {stop("alpha must be > 0 and <= 1")}
  enuzun<-max(length(p),length(alpha))
  p<-rep(p,enuzun/length(p)+1)[1:enuzun]
  alpha<-rep(alpha,enuzun/length(alpha)+1)[1:enuzun]
  qfonk<-NULL
  for (i in 1:enuzun)
  {
    qfonk[i]<-(-pracma::lambertWn(-(exp((-1+log(1-p[i])*alpha[i])/alpha[i])/
                        alpha[i]))*alpha[i]-1+log(1-p[i])*alpha[i])/alpha[i]^2
  }
  if(lower.tail==FALSE)
  {
    qfonk[i]<-(-pracma::lambertWn(-(exp((-1+log(1-(1-p[i]))*alpha[i])/alpha[i])/
                      alpha[i]))*alpha[i]-1+log(1-(1-p[i]))*alpha[i])/alpha[i]^2
  }
  return(qfonk)
}
#' Muth Distribution
#' @export
#' @rdname omd
#' @examples
#' romd(10,alpha=1)
romd<-function(n,alpha)
{
  n<-floor(n)
  if(any(n<1)) {stop("n must be >= 1")}
  if(any(alpha<=0)|any(alpha>1)) {stop("alpha must be > 0 and <= 1")}
  rn<-qomd(stats::runif(n),alpha)
  return(rn)
}
