#' Gamma-Lomax Distribution
#' @export
#' @name gld
#' @param x,q vector of quantiles.
#' @param a,alpha are shape parameters.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#'  to be the number required.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P\left[ X\leq x\right]}, otherwise, \eqn{P\left[ X>x\right] }.
#' @param beta a scale parameter.
#' @description
#' Density, distribution function, quantile function and random generation for
#' the gamma-Lomax distribution with parameters \code{shapes} and \code{scale}.
#' @return \code{dgld} gives the density, \code{pgld} gives the distribution
#' function, \code{qgld} gives the quantile function and \code{rgld} generates
#' random deviates.
#' @details
#' The Gamma-Lomax distribution \code{shape} parameters
#' \eqn{a} and \eqn{\alpha}, and \code{scale} parameter is \eqn{\beta},
#'  has density
#' \deqn{f\left( x\right) =\frac{\alpha \beta ^{\alpha }}
#' {\Gamma \left( a\right)\left( \beta +x\right) ^{\alpha +1}}\left\{ -\alpha
#' \log \left( \frac{\beta }{\beta +x}\right) \right\} ^{a-1},}
#' where
#' \deqn{x>0,~a,\alpha ,\beta >0.}
#' @references  Cordeiro, G. M., Ortega, E. M. and Popović, B. V., 2015,
#' *The gamma-Lomax distribution*, Journal of statistical computation and
#' simulation, 85 (2), 305-319.
#'
#' Ristić, M. M., and Balakrishnan, N. (2012), *The gamma-exponentiated*
#' *exponential distribution. Journal of statistical computation and simulation*
#' , 82(8), 1191-1206.
#' @examples
#' library(new.dist)
#' dgld(1, a=2, alpha=3, beta=4)
dgld<-function(x,a,alpha,beta=1,log=FALSE)
{
  if(any(a<=0)) {stop("a must be > 0")}
  if(any(alpha<=0)) {stop("alpha must be > 0")}
  if(any(beta<=0)) {stop("beta must be > 0")}
  enuzun <- max(length(x),length(a),length(alpha),length(beta))
  x<-rep(x,enuzun/length(x)+1)[1:enuzun]
  a<-rep(a,enuzun/length(a)+1)[1:enuzun]
  alpha<-rep(alpha, enuzun/length(alpha)+1)[1:enuzun]
  beta<-rep(beta,enuzun/length(beta)+1)[1:enuzun]
  pdf<-NULL
  for (i in 1:enuzun)
  {
    suppressWarnings({
    if(x[i]<=0) {pdf[i]<-0} else
    {pdf[i]<-(alpha[i]*beta[i]^alpha[i])/(gamma(a[i])*
        (beta[i]+x[i])^(alpha[i]+1))*(-alpha[i]*log(beta[i]/
                                                  (beta[i]+x[i])))^(a[i]-1)}})
  }
  if(log==TRUE) pdf<-log(pdf)
  return(pdf)
}
#' Gamma-Lomax Distribution
#' @export
#' @rdname gld
#' @examples
#' pgld(1, a=2,alpha=3,beta=4)
pgld<-function(q,a,alpha,beta=1,lower.tail=TRUE,log.p=FALSE)
{
  if(any(a<=0)) {stop("a must be > 0")}
  if(any(alpha<=0)) {stop("alpha must be > 0")}
  if(any(beta<=0)) {stop("beta must be > 0")}
  enuzun <- max(length(q),length(a),length(alpha),length(beta))
  q<-rep(q,enuzun/length(q)+1)[1:enuzun]
  a<-rep(a,enuzun/length(a)+1)[1:enuzun]
  alpha<-rep(alpha, enuzun/length(alpha)+1)[1:enuzun]
  beta<-rep(beta,enuzun/length(beta)+1)[1:enuzun]
  cdf<-NULL
  for (i in 1:enuzun)
  {
    if (q[i]>0) cdf[i]<-(gamma(a[i])-(expint::gammainc(a[i],(-alpha[i]*
                    log(beta[i]/(beta[i]+q[i]))))))/gamma(a[i]) else cdf[i]<-0
  }
  if(lower.tail==FALSE) cdf<-1-cdf
  if(log.p==TRUE) cdf<-log(cdf)
  return(cdf)
}
#' Gamma-Lomax Distribution
#' @export
#' @rdname gld
#' @examples
#' qgld(.8, a=2,alpha=3,beta=4)
qgld<-function(p,a,alpha,beta=1,lower.tail=TRUE)
{
  if(any(p<0)|any(p>1)) {stop("p must be between >= 0 and <= 1")}
  if(any(a<=0)) {stop("a must be > 0")}
  if(any(alpha<=0)) {stop("alpha must be > 0")}
  if(any(beta<=0)) {stop("beta must be > 0")}
  enuzun <- max(length(p),length(a),length(alpha),length(beta))
  p<-rep(p,enuzun/length(p)+1)[1:enuzun]
  a<-rep(a,enuzun/length(a)+1)[1:enuzun]
  alpha<-rep(alpha,enuzun/length(alpha)+1)[1:enuzun]
  beta<-rep(beta,enuzun/length(beta)+1)[1:enuzun]
  kok<-NULL
  quant<-NULL
  for (i in 1:enuzun){
  quant[i]<-VGAM::qlomax(1-exp(-(stats::qgamma(p[i],a[i],1))),beta[i],alpha[i])}
  return(quant)
}
#' Gamma-Lomax Distribution
#' @export
#' @rdname gld
#' @examples
#' rgld(10, a=2,alpha=3,beta=4)
rgld<-function(n,a,alpha,beta=1)
  {
    n<-floor(n)
    if(any(n<1)) {stop("n must be >= 1")}
    if(any(a<=0)) {stop("a must be > 0")}
    if(any(alpha<=0)) {stop("alpha must be > 0")}
    if(any(beta<=0)) {stop("beta must be > 0")}
    suppressWarnings({
    rn<-qgld(stats::runif(n),a,alpha,beta)})
    return(rn)
  }
