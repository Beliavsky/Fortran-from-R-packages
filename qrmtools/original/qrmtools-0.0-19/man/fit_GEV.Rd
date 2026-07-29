\name{fit_GEV}
\alias{fit_GEV_quantile}
\alias{fit_GEV_PWM}
\alias{logLik_GEV}
\alias{fit_GEV_MLE}
\title{Parameter Estimators of the Generalized Extreme Value Distribution}
\description{
  Quantile matching estimator, probability weighted moments estimator,
  log-likelihood and maximum-likelihood estimator for the parameters of
  the generalized extreme value distribution (GEV).
}
\usage{
fit_GEV_quantile(x, p = c(0.25, 0.5, 0.75), cutoff = 3)
fit_GEV_PWM(x)

logLik_GEV(param, x)
fit_GEV_MLE(x, init = c("shape0", "PWM", "quantile"),
            estimate.cov = TRUE, control = list(), ...)
}
\arguments{
  \item{x}{numeric vector of data. In the block maxima method, these are the
    block maxima.}
  \item{p}{\code{numeric(3)} specifying the probabilities whose quantiles are
    matched.}
  \item{cutoff}{positive \eqn{z} after which \eqn{\exp(-z)}{exp(-z)}
    is truncated to 0.}
  \item{param}{\code{numeric(3)} containing the value of the shape
    \eqn{\xi}{xi} (a real), location \eqn{\mu}{mu} (a real) and scale
    \eqn{\sigma}{sigma} (positive real) parameters of the GEV
    distribution in this order.}
  \item{init}{\code{character} string specifying the method for
    computing initial values. Can also be \code{numeric(3)}
    for directly providing \eqn{\xi}{xi}, \eqn{\mu}{mu}, \eqn{\sigma}{sigma}.}
  \item{estimate.cov}{\code{\link{logical}} indicating whether the
    asymptotic covariance matrix of the parameter estimators is to be
    estimated (inverse of observed Fisher information (negative Hessian
    of log-likelihood evaluated at MLE)) and standard errors for the
    estimators of \eqn{\xi}{xi}, \eqn{\mu}{mu}, \eqn{\sigma}{sigma}
    returned, too.}
  \item{control}{\code{\link{list}}; passed to the underlying
    \code{\link{optim}()}.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{optim}()}.}
}
\value{
  \code{fit_GEV_quantile()} and \code{fit_GEV_PWM()} return a
  \code{numeric(3)} giving the parameter estimates for the GEV
  distribution.

  \code{logLik_GEV()} computes the log-likelihood of the GEV
  distribution (\code{-Inf} if not admissible).

  \code{fit_GEV_MLE()} returns the return object of \code{\link{optim}()}
  (by default, the return value \code{value} is the log-likelihood) and,
  appended, the estimated asymptotic covariance matrix and
  standard errors of the parameter estimators, if \code{estimate.cov}.
}
\details{
  \code{fit_GEV_quantile()} matches the empirical \code{p}-quantiles.

  \code{fit_GEV_PWM()} computes the probability weighted moments (PWM)
  estimator of Hosking et al. (1985); see also Landwehr and Wallis (1979).

  \code{fit_GEV_MLE()} uses, as default, the case \eqn{\xi = 0}{xi = 0}
  for computing initial values; this is actually a small positive value
  since Nelder--Mead could fail otherwise. For the other available
  methods for computing initial values, \eqn{\sigma}{sigma} (obtained
  from the case \eqn{\xi = 0}{xi = 0}) is doubled in order to guarantee
  a finite log-likelihood at the initial values. After several
  experiments (see the source code), one can safely say that finding
  initial values for fitting GEVs via MLE is non-trivial; see also the
  block maxima method script about the Black Monday event on
  \url{https://qrmtutorial.org}.

  Caution: See Coles (2001, p. 55) for how to interpret \eqn{\xi\le
    -0.5}{xi <= -0.5}; in particular, the standard asymptotic properties
  of the MLE do not apply.
}
\author{Marius Hofert}
\references{
   McNeil, A. J., Frey, R. and Embrechts, P. (2015).
   \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
   Princeton University Press.

   Hosking, J. R. M., Wallis, J. R. and Wood, E. F. (1985).
   Estimation of the Generalized Extreme-Value Distribution by the
   Method of Probability-Weighted Moments. \emph{Technometrics}
   \bold{27}(3), 251--261.

   Landwehr, J. M. and Wallis, J. R. (1979).
   Probability Weighted Moments Compared With Some Traditional
   Techniques in Estimating Gumbel Parameters and Quantiles.
   \emph{Water Resourches Research} \bold{15}(5), 1055--1064.

   Coles, S. (2001).
   \emph{An Introduction to Statistical Modeling of Extreme Values}.
   Springer-Verlag.
}
\examples{
## Simulate some data
xi <- 0.5
mu <- -2
sig <- 3
n <- 1000
set.seed(271)
X <- rGEV(n, shape = xi, loc = mu, scale = sig)

## Fitting via matching quantiles
(fit.q <- fit_GEV_quantile(X))
stopifnot(all.equal(fit.q[["shape"]], xi,  tol = 0.12),
          all.equal(fit.q[["loc"]],   mu,  tol = 0.12),
          all.equal(fit.q[["scale"]], sig, tol = 0.005))

## Fitting via PWMs
(fit.PWM <- fit_GEV_PWM(X))
stopifnot(all.equal(fit.PWM[["shape"]], xi,  tol = 0.16),
          all.equal(fit.PWM[["loc"]],   mu,  tol = 0.15),
          all.equal(fit.PWM[["scale"]], sig, tol = 0.08))

## Fitting via MLE
(fit.MLE <- fit_GEV_MLE(X))
(est <- fit.MLE$par) # estimates of xi, mu, sigma
stopifnot(all.equal(est[["shape"]], xi,  tol = 0.07),
          all.equal(est[["loc"]],   mu,  tol = 0.12),
          all.equal(est[["scale"]], sig, tol = 0.06))
fit.MLE$SE # estimated asymp. variances of MLEs = std. errors of MLEs

## Plot the log-likelihood in the shape parameter xi for fixed
## location mu and scale sigma (fixed as generated)
xi. <- seq(-0.1, 0.8, length.out = 65)
logLik <- sapply(xi., function(xi..) logLik_GEV(c(xi.., mu, sig), x = X))
plot(xi., logLik, type = "l", xlab = expression(xi),
     ylab = expression("GEV distribution log-likelihood for fixed"~mu~"and"~sigma))
## => Numerically quite challenging (for this seed!)

## Plot the profile likelihood for these xi's
## Note: As initial values for the nuisance parameters mu, sigma, we
##       use their values in the case xi = 0 (for all fixed xi = xi.,
##       in particular those xi != 0). Furthermore, for the given data X
##       and xi = xi., we make sure the initial value for sigma is so large
##       that the density is not 0 and thus the log-likelihood is finite.
pLL <- sapply(xi., function(xi..) {
    scale.init <- sqrt(6 * var(X)) / pi
    loc.init <- mean(X) - scale.init * 0.5772157
    while(!is.finite(logLik_GEV(c(xi.., loc.init, scale.init), x = X)) &&
          is.finite(scale.init)) scale.init <- scale.init * 2
    optim(c(loc.init, scale.init), fn = function(nuis)
                logLik_GEV(c(xi.., nuis), x = X),
    		        control = list(fnscale = -1))$value
})
plot(xi., pLL, type = "l", xlab = expression(xi),
     ylab = "GEV distribution profile log-likelihood")
}
\keyword{distribution}
