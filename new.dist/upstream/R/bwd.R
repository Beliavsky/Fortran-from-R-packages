#' Bimodal Weibull Distribution
#' @export
#' @name bwd
#' @param x,q vector of quantiles.
#' @param alpha a shape parameter.
#' @param beta a scale parameter.
#' @param sigma a control parameter that controls the uni- or bimodality of the
#' distribution.
#' @param p vector of probabilities.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken
#' to be the number required.
#' @param log,log.p logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail logical; if TRUE (default), probabilities are
#' \eqn{P\left[ X\leq x\right]}, otherwise, \eqn{P\left[ X>x\right] }.
#' @description
#' Density, distribution function, quantile function and random generation for
#' a Bimodal Weibull distribution with parameters \code{shape} and \code{scale}.
#' @return \code{dbwd} gives the density, \code{pbwd} gives the distribution
#' function, \code{qbwd} gives the quantile function and \code{rbwd} generates
#' random deviates.
#' @details
#' A Bimodal Weibull distribution with \code{shape} parameter \eqn{\alpha},
#' \code{scale} parameter \eqn{\beta},and the \code{control} parameter
#' \eqn{\sigma} that determines the uni- or bimodality of the
#' distribution, has density
#' \deqn{f\left( x\right) =\frac{\alpha }{\beta Z_{\theta }}
#' \left[ 1+\left( 1-\sigma~x\right) ^{2}\right] \left( \frac{x}{\beta }
#' \right) ^{\alpha -1}\exp \left( -\left( \frac{x}{\beta }\right) ^{\alpha }
#' \right),}
#' where
#' \deqn{Z_{\theta }=2+\sigma ^{2}\beta ^{2}\Gamma
#' \left( 1+\left( 2/\alpha \right)\right) -2\sigma \beta \Gamma
#' \left( 1+\left( 1/\alpha \right) \right) }
#' and
#' \deqn{x\geq 0,~\alpha ,\beta >0~ and ~\sigma \in\mathbb{R}.}
#' @references  Vila, R. and Niyazi Çankaya, M., 2022,
#' *A bimodal Weibull distribution: properties and inference*,
#' Journal of Applied Statistics, 49 (12), 3044-3062.
#' @examples
#' library(new.dist)
#' dbwd(1,alpha=2:10,beta=3,sigma=4)
dbwd<-function(x,alpha,beta=1,sigma, log = FALSE)
{
  if(any(alpha<=0)) {stop("alpha must be > 0")}
  if(any(beta<=0)) {stop("beta must be > 0")}
  enuzun<-max(length(x),length(alpha),length(beta),length(sigma))
  x<-rep(x,enuzun/length(x)+1)[1:enuzun]
  alpha<-rep(alpha,enuzun/length(alpha)+1)[1:enuzun]
  beta<-rep(beta,enuzun/length(beta)+1)[1:enuzun]
  sigma<-rep(sigma,enuzun/length(sigma)+1)[1:enuzun]
  pdf<-NULL
  z<-NULL
  for (i in 1:enuzun)suppressWarnings(
    {
      if(x[i]<0) {pdf[i]<-0} else
      {z[i]<-2+sigma[i]^2*beta[i]^2*gamma(1+(2/alpha[i]))-2*sigma[i]*beta[i]*
        gamma(1+(1/alpha[i]))
      pdf[i]<-(alpha[i]/(beta[i]*z[i]))*(1+(1-sigma[i]*x[i])^2)*
        (x[i]/beta[i])^(alpha[i]-1)*exp(-(x[i]/beta[i])^alpha[i])
      }
    })
  if(log) {pdf<-log(pdf)}
  return(pdf)
}
#' Bimodal Weibull Distribution
#' @export
#' @rdname bwd
#' @examples
#' pbwd(1,alpha=2,beta=3,sigma=4)
pbwd<-function(q,alpha,beta=1,sigma,lower.tail=TRUE,log.p=FALSE)
{
  if(any(alpha<=0)) {stop("alpha must be > 0")}
  if(any(beta<=0)) {stop("beta must be > 0")}
  enuzun<-max(length(q),length(alpha),length(beta),length(sigma))
  q<-rep(q,enuzun/length(q)+1)[1:enuzun]
  alpha<-rep(alpha,enuzun/length(alpha)+1)[1:enuzun]
  beta<-rep(beta,enuzun/length(beta)+1)[1:enuzun]
  sigma<-rep(sigma,enuzun/length(sigma)+1)[1:enuzun]
  cdf<-NULL
  for (i in 1:enuzun)suppressWarnings(
    {
      if(q[i]>=0) (cdf[i]<-(2-(1+(1-sigma[i]*q[i])^2)*
            exp(-(q[i]/beta[i])^alpha[i])-((2*sigma[i]*beta[i])/alpha[i])*
        (pracma::gammainc(q[i]^alpha[i]/beta[i]^alpha[i],1/alpha[i])[1]-sigma[i]
                        *beta[i]*pracma::gammainc(q[i]^alpha[i]/beta[i]^
        alpha[i],2/alpha[i])[1]))/(2+sigma[i]^2*beta[i]^2*gamma(1+(2/alpha[i]))
                    -2*sigma[i]*beta[i]*gamma(1+(1/alpha[i])))) else cdf[i]<-0
    })
  if(!lower.tail) cdf<-1-cdf
  if(log.p) cdf<-log(cdf)
  return(cdf)
}
#' Bimodal Weibull Distribution
#' @export
#' @rdname bwd
#' @examples
#' qbwd(.7,alpha=2,beta=3,sigma=4)
qbwd<-function(p,alpha,beta=1,sigma,lower.tail=TRUE)
{
  if(any(p<0)|any(p>1)) {stop("p must be between >= 0 and <= 1")}
  if(any(alpha<=0)) {stop("alpha must be > 0")}
  if(any(beta<=0)) {stop("beta must be > 0")}
  enuzun<-max(length(p),length(alpha),length(beta),length(sigma))
  p<-rep(p,enuzun/length(p)+1)[1:enuzun]
  alpha<-rep(alpha,enuzun/length(alpha)+1)[1:enuzun]
  beta<-rep(beta,enuzun/length(beta)+1)[1:enuzun]
  sigma<-rep(sigma,enuzun/length(sigma)+1)[1:enuzun]
  kok<-NULL
  for (i in 1:enuzun)suppressWarnings(
    {
      Ex<-beta[i]*(2*gamma(1+(1/alpha[i]))+sigma[i]^2*beta[i]^2*
                     gamma(1+(3/alpha[i]))-2*sigma[i]*beta[i]*
                     gamma(1+(2/alpha[i])))/
        (2+sigma[i]^2*beta[i]^2*gamma(1+(2/alpha[i]))-2*sigma[i]*beta[i]*
           gamma(1+(1/alpha[i])))
      Y<-function(x)
      {
        abs(((2-(1+(1-sigma[i]*x)^2)*exp(-(x/beta[i])^alpha[i])-
                ((2*sigma[i]*beta[i])/alpha[i])*
                (pracma::gammainc(x^alpha[i]/beta[i]^
                                    alpha[i],1/alpha[i])[1]-sigma[i]*beta[i]*
                  pracma::gammainc(x^alpha[i]/beta[i]^alpha[i],2/alpha[i])[1]))
             /(2+sigma[i]^2*beta[i]^2*gamma(1+(2/alpha[i]))-2*sigma[i]*beta[i]*
                 gamma(1+(1/alpha[i]))))-p[i])
      }
      if(lower.tail==FALSE)
      {
        Y<-function(x)
        {
          abs(((2-(1+(1-sigma[i]*x)^2)*exp(-(x/beta[i])^alpha[i])-((2*sigma[i]*
              beta[i])/alpha[i])*(pracma::gammainc(x^alpha[i]/beta[i]^alpha[i],
                        1/alpha[i])[1]-sigma[i]*beta[i]*
      pracma::gammainc(x^alpha[i]/beta[i]^alpha[i],2/alpha[i])[1]))/
        (2+sigma[i]^2*beta[i]^2*gamma(1+(2/alpha[i]))-2*sigma[i]*beta[i]*
                                              gamma(1+(1/alpha[i]))))-(1-p[i]))
        }
      }
      kok[i]<-stats::optim(Ex,Y)$par
    })
  return(kok)
}
#' Bimodal Weibull Distribution
#' @export
#' @rdname bwd
#' @examples
#' rbwd(10,alpha=2,beta=3,sigma=4)
rbwd<-function(n,alpha,beta=1,sigma)
{
  n<-floor(n)
  if(any(n<1)) {stop("n must be >= 1")}
 if(any(alpha<=0)) {stop("alpha must be > 0")}
  if(any(beta<=0)) {stop("beta must be > 0")}
  suppressWarnings({
    rn<-qbwd(stats::runif(n),alpha,beta,sigma)})
  return(rn)
}
