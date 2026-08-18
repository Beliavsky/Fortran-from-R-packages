#' @title LindleyPoisson
#' @aliases plindleypoisson dlindleypoisson hlindleypoisson qlindleypoisson rlindleypoisson
#' @description distribution function, density function, hazard rate function, quantile function, random number generation
#' @author Saralees Nadarajah & Yuancheng Si \email{siyuanchengman@gmail.com}
#' @author Peihao Wang
#' @param x vector of positive quantiles.
#' @param lambda positive parameter
#' @param theta positive parameter.
#' @param n number of observations.
#' @param p vector of probabilities.
#' @param log.p logical; If \code{TRUE}, probabilities \eqn{p} are given as \eqn{log(p)}.
#' @references
#'
#' Si, Y. & Nadarajah, S., (2018). Lindley Power Series Distributions. \emph{Sankhya A}, \bold{9}, pp1-15.
#'
#' Ghitany, M. E., Atieh, B., Nadarajah, S., (2008). Lindley distribution and its application. \emph{Mathematics and Computers in Simulation}, \bold{78}, (4), 49-506.
#'
#' Jodra, P., (2010). Computer generation of random variables with Lindley or Poisson-Lindley distribution via the Lambert W function. \emph{Mathematics and Computers in Simulation}, \bold{81}, (4), 851-859.
#'
#' Lindley, D. V., (1958). Fiducial distributions and Bayes' theorem. \emph{Journal of the Royal Statistical Society. Series B. Methodological}, \bold{20}, 102-107.
#'
#' Lindley, D. V., (1965). \emph{Introduction to Probability and Statistics from a Bayesian View-point, Part II: Inference}. Cambridge University Press, New York.
#'
#' @details
#'
#' Probability density function
#' \deqn{f(x)=\frac{\theta\lambda^2}{(\lambda+1)A(\theta)}(1+x)exp(-\lambda x)A^{'}(\phi)}
#'
#' Cumulative distribution function
#' \deqn{F(x)=\frac{A(\phi)}{A(\theta)}}
#'
#' Quantile function
#' \deqn{F^{-1}(p)=-1-\frac{1}{\lambda}-\frac{1}{\lambda}W_{-1}\left\{\frac{\lambda+1}{exp(\lambda+1)}\left[\frac{1}{\theta}A^{-1}\{pA(\theta)\}-1\right]\right\}}
#'
#' Hazard rate function
#' \deqn{h(x)=\frac{\theta\lambda^2}{1+\lambda}(1+x)exp(-\lambda x)\frac{A^{'}(\phi)}{A(\theta)-A(\phi)}}
#'
#' where \eqn{W_{-1}} denotes the negative branch of the Lambert W function. \eqn{A(\theta)=\sum_{n=1}^{\infty}a_n\theta^{n}} is given by specific power series distribution.
#' Note that \eqn{x>0, \lambda>0}  for all members in Lindley Power Series distribution.
#' \eqn{0<\theta<1}  for Lindley-Geometric distribution, Lindley-logarithmic distribution, Lindley-Negative Binomial distribution.
#' \eqn{\theta>0} for Lindley-Poisson distribution, Lindley-Binomial distribution.
#' @return \code{plindleypoisson} gives the culmulative distribution function
#' @return \code{dlindleypoisson} gives the probability density function
#' @return \code{hlindleypoisson} gives the hazard rate function
#' @return \code{qlindleypoisson} gives the quantile function
#' @return \code{rlindleypoisson} gives the random number generatedy by distribution
#' @return Invalid arguments will return an error message.
#'
#' @examples
#' set.seed(1)
#' lambda = 1
#' theta = 0.5
#' n = 10
#' x <- seq(from = 0.1,to = 6,by = 0.5)
#' p <- seq(from = 0.1,to = 1,by = 0.1)
#' plindleypoisson(x, lambda, theta, log.p = FALSE)
#' dlindleypoisson(x, lambda, theta)
#' hlindleypoisson(x, lambda, theta)
#' qlindleypoisson(p, lambda, theta)
#' rlindleypoisson(n, lambda, theta)
#' @rdname LindleyPoisson
#' @export
plindleypoisson <- function(x , lambda , theta , log.p = FALSE)
{
  stopifnot(theta > 0,lambda > 0,x > 0,is.logical(log.p))
  t1 = 1 - (lambda + 1 + lambda * x) / (lambda + 1) * exp(-lambda * x)
  cdf = (exp(theta * t1) - 1) / (exp(theta) - 1)
  if(log.p) return(log(cdf)) else return(cdf)
}



#' @rdname LindleyPoisson
#' @export
dlindleypoisson <- function(x, lambda, theta)
{
  stopifnot(theta > 0,lambda > 0,x > 0)
  t0 = theta * (1 - (lambda + 1 + lambda * x) / (1 + lambda) * exp(-lambda * x))
  t1 = theta * lambda ** 2 / ((lambda + 1) * (exp(theta) - 1)) * (1 + x) * exp(-lambda * x)
  pdf = t1 * exp(t0)
  return(pdf)
}


#' @rdname LindleyPoisson
#' @export
hlindleypoisson <- function(x, lambda, theta)
{
  stopifnot(theta > 0,lambda > 0,x > 0)
  t0 = theta * (1 - (lambda + 1 + lambda * x) / (1 + lambda) * exp(-lambda * x))
  t1 = theta * lambda ** 2 / (lambda + 1) *  (1 + x) * exp(-lambda * x)
  t2 = exp(t0) / (exp(theta) - exp(t0))
  hazard = t0 * t1 * t2
  return(hazard)
}

#' @rdname LindleyPoisson
#' @export
qlindleypoisson <- function(p, lambda, theta)
{
  stopifnot(theta > 0,lambda > 0)
  t0 = p * (exp(theta) - 1)
  t1 = log(t0 + 1) / theta - 1
  t2 = (lambda + 1) / exp(lambda + 1) * t1
  x = - lamW::lambertWm1(t2) / lambda - 1 / lambda -1
  return(x)
}


#' method #1
#' #' @rdname LindleyPoisson
#' #' @export
#' rlindleypoisson <- function(n, lambda, theta)
#' {
#'   stopifnot(theta > 0,lambda > 0,n %% 1==0)
#'   i=0
#'   randdata=c()
#'   while (i<n)
#'   {
#'     N = rpois(1, theta)
#'     data = LindleyR::rlindley(N, lambda, mixture = TRUE)
#'     randdata=c(randdata,max(data))
#'     i=i+1
#'   }
#'   return(randdata)
#'
#'
#' }
#' method #2
#' @rdname LindleyPoisson
#' @export
rlindleypoisson <- function(n, lambda, theta)
{
  stopifnot(theta > 0,lambda > 0,n %% 1==0)
  y=stats::runif(n, min=0, max = 1)
  randdata=qlindleypoisson(y, lambda, theta)
  return(randdata)


}
