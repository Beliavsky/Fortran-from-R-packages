#' Power Log Dagum Distribution
#' @export
#' @name pldd
#' @param x,q vector of quantiles.
#' @param alpha,beta,theta are parameters.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#'  to be the number required.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P\left[ X\leq x\right]}, otherwise, \eqn{P\left[ X>x\right] }.
#' @description
#' Density, distribution function, quantile function and random generation for
#' a Power Log Dagum distribution.
#' @return \code{dpldd} gives the density, \code{ppldd} gives the distribution
#' function, \code{qpldd} gives the quantile function and \code{rpldd} generates
#'  random deviates.
#' @details
#' A Power Log Dagum Distribution with parameters \eqn{\alpha}, \eqn{\beta} and
#' \eqn{\theta}, has density
#'  \deqn{f\left( x\right) =\alpha
#'  \left( \beta +\theta \left\vert x\right\vert^{\beta -1}
#'  \right) e^{-\left( \beta x+sign\left( x\right)
#'  \left( \theta/\beta \right) \left\vert
#'    x\right\vert ^{\beta }\right) ~}~\left(1+e^{-\left( \beta x+sign
#'    \left( x\right)\left( \theta /\beta \right)
#'    \left\vert x\right\vert ^{\beta }\right) }
#'    \right) ^{-\left( \alpha +1\right)},}
#' where
#' \deqn{x\in \mathbb{R},~\beta \in \mathbb{R},~\alpha >0~and~\theta \geq 0}
#' @note
#' The distributions hazard function
#' \deqn{h\left( x\right) =\frac{\alpha
#' \left( \beta +\theta \left\vert x\right\vert^{\beta -1}
#' \right) e^{-\left( \beta x+sign\left( x\right) \left( \theta/\beta \right)
#' \left\vert x\right\vert ^{\beta }\right) }\left( 1+e^{-\left(\beta x+sign
#' \left( x\right) \left( \theta /\beta \right) \left\vert x
#' \right\vert ^{\beta }\right) }\right) ^{-\left(\alpha +1\right) }}
#'         {1-\left( 1+e^{-\left( \beta x+sign\left( x\right) \left( \theta /
#' \beta \right) \left\vert x\right\vert ^{\beta }\right) }
#' \right) ^{-\alpha }} .}
#' @references  Bakouch, H. S., Khan, M. N., Hussain, T. and Chesneau, C., 2019,
#' *A power log-Dagum distribution: estimation and applications*, Journal of
#' Applied Statistics, 46 (5), 874-892.
#' @examples
#' library(new.dist)
#' dpldd(1, alpha=2, beta=3, theta=4)
dpldd<-function(x,alpha,beta,theta,log=FALSE)
{
    if(any(alpha<=0)) {stop("alpha must be > 0")}
    if(any(theta<0)) {stop("theta must be >= 0")}
    enuzun<-max(length(x),length(alpha),length(beta),length(theta))
    x<-rep(x,enuzun/length(x)+1)[1:enuzun]
    alpha<-rep(alpha,enuzun/length(alpha)+1)[1:enuzun]
    beta<-rep(beta,enuzun/length(beta)+1)[1:enuzun]
    theta<-rep(theta,enuzun/length(theta)+1)[1:enuzun]
    pdf<-NULL
    for (i in 1:enuzun)
    {suppressWarnings({
      pdf[i]<-alpha[i]*(beta[i]+theta[i]*abs(x[i])^(beta[i]-1))*(exp(-(beta[i]*
                    x[i]+(sign(x[i]))*(theta[i]/beta[i])*abs(x[i])^beta[i])))*
      (1+exp(-(beta[i]*x[i]+(sign(x[i]))*(theta[i]/beta[i])*
                 abs(x[i])^beta[i])))^(-(alpha[i]+1))})
    }
    if(log==TRUE) pdf<-log(pdf)
    return(pdf)
  }
#' Power Log Dagum Distribution
#' @export
#' @rdname pldd
#' @examples
#' ppldd(1,alpha=2,beta=3,theta=4)
ppldd<-function(q,alpha,beta,theta,lower.tail=TRUE,log.p=FALSE)
  {
  if(any(alpha<=0)) {stop("alpha must be > 0")}
  if(any(theta<0)) {stop("theta must be >= 0")}
    enuzun<-max(length(q),length(alpha),length(beta),length(theta))
    q<-rep(q,enuzun/length(q)+1)[1:enuzun]
    alpha<-rep(alpha,enuzun/length(alpha)+1)[1:enuzun]
    beta<-rep(beta,enuzun/length(beta)+1)[1:enuzun]
    theta<-rep(theta,enuzun/length(theta)+1)[1:enuzun]
    cdf<-NULL
    for (i in 1:enuzun) suppressWarnings(
    {
      cdf[i]<-(1+exp(-(beta[i]*q[i]+(sign(q[i]))*(theta[i]/beta[i])*
                         abs(q[i])^beta[i])))^(-alpha[i])
    })
    if(lower.tail==FALSE) cdf<-1-cdf
    if(log.p==TRUE) cdf<-log(cdf)
    return(cdf)
}
#' Power Log Dagum Distribution
#' @export
#' @rdname pldd
#' @examples
#' qpldd(.8,alpha=2,beta=3,theta=4)
qpldd<-function(p,alpha,beta,theta,lower.tail=TRUE)
  {
  if(any(p<0)|any(p>1)) {stop("p must be between >= 0 and <= 1")}
  if(any(alpha<=0)) {stop("alpha must be > 0")}
  if(any(theta<0)) {stop("theta must be >= 0")}
    enuzun<-max(length(p),length(alpha),length(beta),length(theta))
    p<-rep(p,enuzun/length(p)+1)[1:enuzun]
    alpha<-rep(alpha,enuzun/length(alpha)+1)[1:enuzun]
    beta<-rep(beta,enuzun/length(beta)+1)[1:enuzun]
    theta<-rep(theta,enuzun/length(theta)+1)[1:enuzun]
    kok<-NULL
    for (i in 1:enuzun) suppressWarnings(
    {
      Y<-function(x)
      {
        (1+exp(-(beta[i]*x+(sign(x))*(theta[i]/beta[i])*
                   abs(x)^beta[i])))^(-alpha[i])-p[i]
      }
      if(lower.tail==FALSE)
      {
        Y<-function(x)
        {
          (1+exp(-(beta[i]*x+(sign(x))*(theta[i]/beta[i])*
                     abs(x)^beta[i])))^(-alpha[i])-(1-p[i])
        }
      }
      kok[i]<-(stats::uniroot(Y,c(-10,10)))$root
    })
    return(kok)
  }
#' Power Log Dagum Distribution
#' @export
#' @rdname pldd
#' @examples
#' rpldd(10,alpha=2,beta=3,theta=4)
rpldd<-function(n,alpha,beta,theta)
  {
    n<-floor(n)
    if(any(n<1)) {stop("n must be >= 1")}
    if(any(alpha<=0)) {stop("alpha must be > 0")}
    if(any(theta<0)) {stop("theta must be >= 0")}
    suppressWarnings({
    rn<-qpldd(stats::runif(n),alpha,beta,theta)})
    return(rn)
  }
