#' Maxwell Distribution
#' @export
#' @name md
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
#' Maxwell distribution with parameter \code{scale}.
#' @return \code{dmd} gives the density, \code{pmd} gives the distribution
#' function, \code{qmd} gives the quantile function and \code{rmd} generates
#' random deviates.
#' @details
#' Maxwell distribution with \code{scale} parameter \eqn{\theta},
#' has density
#' \deqn{f\left( x\right) =\frac{4}{\sqrt{\pi }}
#' \frac{1}{\theta ^{3/2}}x^{2}e^{-x^{2}/\theta },}
#' where
#' \deqn{0\leq x<\infty ,~~\theta >0.}
#' @references  Krishna, H., Vivekanand and Kumar, K., 2015,
#' *Estimation in Maxwell distribution with randomly censored data*, Journal of
#' statistical computation and simulation, 85 (17), 3560-3578.
#' @examples
#' library(new.dist)
#' dmd(1,theta=2)
dmd<-function(x,theta=1,log=FALSE)
{
  if(any(theta<=0)) {stop("theta must be > 0")}
  enuzun<-max(length(x),length(theta))
  x<-rep(x,enuzun/length(x)+1)[1:enuzun]
  theta<-rep(theta,enuzun/length(theta)+1)[1:enuzun]
  pdf<-NULL
  for (i in 1:enuzun)
  {
    if(x[i]<0) {pdf[i]<-0} else
    {pdf[i]<-(4/pi^(1/2))*(1/theta[i]^(3/2))*x[i]^2*exp((-x[i]^2)/theta[i])}
  }
  if(log==TRUE) pdf<-log(pdf)
  return(pdf)
}
#' @export
#' @rdname md
#' @examples
#' pmd(1,theta=2)
pmd<-function(q,theta=1,lower.tail=TRUE,log.p=FALSE)
{
  if(any(theta<=0)) {stop("theta must be > 0")}
  enuzun<-max(length(q),length(theta))
  q<-rep(q,enuzun/length(q)+1)[1:enuzun]
  theta<-rep(theta,enuzun/length(theta)+1)[1:enuzun]
  cdf<-NULL
  for (i in 1:enuzun)
  {
    if(q[i]>=0)cdf[i]<-pracma::gammainc((q[i]^2/theta[i]),3/2)[3] else cdf[i]<-0
  }
  if(lower.tail==FALSE) cdf<-1-cdf
  if(log.p==TRUE) cdf<-log(cdf)
  return(cdf)
}
#' @export
#' @rdname md
#' @examples
#' qmd(.4,theta=5)
qmd<-function(p,theta=1,lower.tail=TRUE)
  {
  if(any(p<0)|any(p>1)) {stop("p must be between >= 0 and <= 1")}
  if(any(theta<=0)) {stop("theta must be > 0")}
  enuzun<-max(length(p),length(theta))
    p<-rep(p,enuzun/length(p)+1)[1:enuzun]
    theta<-rep(theta,enuzun/length(theta)+1)[1:enuzun]
    kok<-NULL
    for (i in 1:enuzun)suppressWarnings(
    {
      Ex<-2*(theta[i]/pi)^(1/2)
      Y<-function(x)
      {
        abs(pracma::gammainc((x^2/theta[i]),3/2)[3]-p[i])
      }
      if(lower.tail==FALSE)
      {
        Y<-function(x)
        {
          abs(pracma::gammainc((x^2/theta[i]),3/2)[3]-(1-p[i]))
        }
      }
      kok[i]<-stats::optim(Ex,Y)$par
    })
    return(kok)
  }
#' @export
#' @rdname md
#' @examples
#' rmd(10,theta=1)
rmd<-function(n,theta=1)
  {
    n<-floor(n)
    if(any(n<1)) {stop("n must be >= 1")}
    if(any(theta<=0)) {stop("theta must be > 0")}
    suppressWarnings({
    rn<-qmd(stats::runif(n),theta)
    return(rn)})
  }

