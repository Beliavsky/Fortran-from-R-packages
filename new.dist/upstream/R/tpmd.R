#' Power Muth Distribution
#' @export
#' @name tpmd
#' @param x,q vector of quantiles.
#' @param beta a scale parameter.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#'  to be the number required.
#' @param alpha a shape parameter.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P\left[ X\leq x\right]}, otherwise, \eqn{P\left[ X>x\right] }.
#' @description
#' Density, distribution function, quantile function and random generation for
#' the Power Muth distribution with parameters \code{shape} and \code{scale}.
#' @return \code{dtpmd} gives the density, \code{ptpmd} gives the distribution
#' function, \code{qtpmd} gives the quantile function and \code{rtpmd} generates
#'  random deviates.
#' @details
#' The Power Muth distribution with \code{shape} parameter \eqn{\alpha} and
#' \code{scale} parameter \eqn{\beta} has density
#' \deqn{f\left( x\right) =\frac{\alpha }{\beta ^\alpha }x^{\alpha -1}
#' \left( e^{\left(x/\beta \right) ^{\alpha }}-1\right)
#' \left( e^{\left( x/\beta \right)  ^{\alpha }-
#' \left( e^{\left( x/\beta \right) ^{\alpha }}-1\right) }\right), }
#' where
#' \deqn{x>0,~\alpha ,\beta>0.}
#' @references  Jodra, P., Gomez, H. W., Jimenez-Gamero,
#'  M. D., and Alba-Fernandez, M. V. (2017).
#' *The power Muth distribution* . Mathematical Modelling and Analysis, 22(2),
#'  186-201.
#' @note
#' Hazard function;
#' \deqn{h\left( \beta ,\alpha \right) =\frac{\alpha }{\beta ^{\alpha }}
#' \left(e^{\left( x/\beta \right) ^{\alpha }}-1\right) x^{\alpha -1}}
#' @examples
#' library(new.dist)
#' dtpmd(1, beta=2, alpha=3)
dtpmd<-function(x,beta=1,alpha,log=FALSE)
{
  if(any(beta<=0)) {stop("beta must be > 0")}
  if(any(alpha<=0)) {stop("alpha must be > 0")}
  enuzun <- max(length(x), length(alpha), length(beta))
  x<-rep(x,enuzun/length(x)+1)[1:enuzun]
  alpha<-rep(alpha, enuzun/length(alpha)+1)[1:enuzun]
  beta<-rep(beta,enuzun/length(beta)+1)[1:enuzun]
  pdf<-NULL
  for (i in 1:enuzun)
  {
    if(x[i]<=0) {pdf[i]<-0} else
    {
      pdf[i]<-(alpha[i]/(beta[i]^alpha[i]))*(x[i]^(alpha[i]-1))*
        (exp((x[i]/beta[i])^alpha[i])-1)*exp(((x[i]/beta[i])^alpha[i])-
                                               (exp((x[i]/beta[i])^alpha[i])-1))
    }
  }
  if(log==TRUE) pdf<-log(pdf)
  return(pdf)
}
#' Power Muth Distribution
#' @export
#' @rdname tpmd
#' @examples
#' ptpmd(1,beta=2,alpha=3)
ptpmd<-function(q,beta=1,alpha,lower.tail=TRUE,log.p=FALSE)
{
  if(any(beta<=0)) {stop("beta must be > 0")}
  if(any(alpha<=0)) {stop("alpha must be > 0")}
  enuzun<-max(length(q), length(alpha), length(beta))
  q<-rep(q,enuzun/length(q)+1)[1:enuzun]
  alpha<-rep(alpha, enuzun/length(alpha)+1)[1:enuzun]
  beta<-rep(beta,enuzun/length(beta)+1)[1:enuzun]
  cdf<-NULL
  for (i in 1:enuzun)
  {
    if(q[i]>0) cdf[i]<-1-exp(((q[i]/beta[i])^alpha[i])-(exp((q[i]/
                                        beta[i])^alpha[i])-1)) else cdf[i]<-0
  }
  if(lower.tail==FALSE) cdf<-1-cdf
  if(log.p==TRUE) cdf<-log(cdf)
  return(cdf)
}
#' Power Muth Distribution
#' @export
#' @rdname tpmd
#' @examples
#' qtpmd(.5,beta=2,alpha=3)
qtpmd<-function(p,beta=1,alpha,lower.tail=TRUE)
{
  if(any(p<0)|any(p>1)) {stop("p must be between >= 0 and <= 1")}
  if(any(beta<=0)) {stop("beta must be > 0")}
  if(any(alpha<=0)) {stop("alpha must be > 0")}
  enuzun <- max(length(p), length(alpha), length(beta))
  p<-rep(p,enuzun/length(p)+1)[1:enuzun]
  alpha<-rep(alpha, enuzun/length(alpha)+1)[1:enuzun]
  beta<-rep(beta,enuzun/length(beta)+1)[1:enuzun]
  quant<-NULL
  for (i in 1:enuzun)
  {
    quant[i]<-exp(log(-pracma::lambertWn((-1+p[i])/exp(1))-1+log(1-p[i]))/
                    alpha[i])*beta[i]
  }
  if(lower.tail==FALSE)
  {
    quant[i]<-(exp(log(-pracma::lambertWn((-1+(1-p[i]))/exp(1))-1+
                         log(1-(1-p[i])))/alpha[i])*beta[i])
  }
  {
    return(quant)
  }
}
#' Power Muth Distribution
#' @export
#' @rdname tpmd
#' @examples
#' rtpmd(10,beta=2,alpha=3)
rtpmd<-function(n,beta=1,alpha)
  {
    n<-floor(n)
    if(any(n<1)) {stop("n must be >= 1")}
    if(any(beta<=0)) {stop("beta must be > 0")}
    if(any(alpha<=0)) {stop("alpha must be > 0")}
    suppressWarnings({
    rn<-qtpmd(stats::runif(n),beta,alpha)})
    return(rn)
  }
