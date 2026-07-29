\name{fit_GPD}
\alias{fit_GPD_MOM}
\alias{fit_GPD_PWM}
\alias{logLik_GPD}
\alias{fit_GPD_MLE}
\title{Parameter Estimators of the Generalized Pareto Distribution}
\description{
  Method-of-moments estimator, probability weighted moments estimator,
  log-likelihood and maximum-likelihood estimator for the parameters of
  the generalized Pareto distribution (GPD).
}
\usage{
fit_GPD_MOM(x)
fit_GPD_PWM(x)

logLik_GPD(param, x)
fit_GPD_MLE(x, init = c("PWM", "MOM", "shape0"),
            estimate.cov = TRUE, control = list(), ...)
}
\arguments{
  \item{x}{numeric vector of data. In the peaks-over-threshold method,
    these are the excesses (exceedances minus threshold).}
  \item{param}{\code{numeric(2)} containing the value of the
    shape \eqn{\xi}{xi} (a real) and scale \eqn{\beta}{beta}
    (positive real) parameters of the GPD in this order.}
  \item{init}{\code{character} string specifying the method for
    computing initial values. Can also be \code{numeric(2)}
    for directly providing \eqn{\xi}{xi} and \eqn{\beta}{beta}.}
  \item{estimate.cov}{\code{\link{logical}} indicating whether the
    asymptotic covariance matrix of the parameter estimators is to be
    estimated (inverse of observed Fisher information (negative Hessian
    of log-likelihood evaluated at MLE)) and standard errors for the
    estimators of \eqn{\xi}{xi} and \eqn{\beta}{beta}
    returned, too.}
  \item{control}{\code{\link{list}}; passed to the underlying
    \code{\link{optim}()}.}
  \item{\dots}{additional arguments passed to the underlying
    \code{\link{optim}()}.}
}
\value{
  \code{fit_GEV_MOM()} and \code{fit_GEV_PWM()} return a
  \code{numeric(3)} giving the parameter estimates for the GPD.

  \code{logLik_GPD()} computes the log-likelihood of the GPD
  (\code{-Inf} if not admissible).

  \code{fit_GPD_MLE()} returns the return object of \code{\link{optim}()}
  and, appended, the estimated asymptotic covariance matrix and
  standard errors of the parameter estimators, if \code{estimate.cov}.
}
\details{
  \code{fit_GPD_MOM()} computes the method-of-moments (MOM) estimator.

  \code{fit_GPD_PWM()} computes the probability weighted moments (PWM)
  estimator of Hosking and Wallis (1987); see also Landwehr et al. (1979).

  \code{fit_GPD_MLE()} uses, as default, \code{fit_GPD_PWM()} for
  computing initial values. The former requires the data \code{x}
  to be non-negative and adjusts \eqn{\beta}{beta} if \eqn{\xi}{xi} is
  negative, so that the log-likelihood at the initial value should be
  finite.
}
\author{Marius Hofert}
\references{
   McNeil, A. J., Frey, R. and Embrechts, P. (2015).
   \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
   Princeton University Press.

   Hosking, J. R. M. and Wallis, J. R. (1987).
   Parameter and Quantile Estimation for the Generalized Pareto Distribution.
   \emph{Technometrics} \bold{29}(3), 339--349.

   Landwehr, J. M., Matalas, N. C. and Wallis, J. R. (1979).
   Estimation of Parameters and Quantiles of Wakeby Distributions.
   \emph{Water Resourches Research} \bold{15}(6), 1361--1379.
}
\examples{
## Simulate some data
xi <- 0.5
beta <- 3
n <- 1000
set.seed(271)
X <- rGPD(n, shape = xi, scale = beta)

## Fitting via matching moments
(fit.MOM <- fit_GPD_MOM(X))
stopifnot(all.equal(fit.MOM[["shape"]], xi,   tol = 0.52),
          all.equal(fit.MOM[["scale"]], beta, tol = 0.24))

## Fitting via PWMs
(fit.PWM <- fit_GPD_PWM(X))
stopifnot(all.equal(fit.PWM[["shape"]], xi,   tol = 0.2),
          all.equal(fit.PWM[["scale"]], beta, tol = 0.12))

## Fitting via MLE
(fit.MLE <- fit_GPD_MLE(X))
(est <- fit.MLE$par) # estimates of xi, mu, sigma
stopifnot(all.equal(est[["shape"]], xi,   tol = 0.12),
          all.equal(est[["scale"]], beta, tol = 0.11))
fit.MLE$SE # estimated asymp. variances of MLEs = std. errors of MLEs

## Plot the log-likelihood in the shape parameter xi for fixed
## scale beta (fixed as generated)
xi. <- seq(-0.1, 0.8, length.out = 65)
logLik <- sapply(xi., function(xi..) logLik_GPD(c(xi.., beta), x = X))
plot(xi., logLik, type = "l", xlab = expression(xi),
     ylab = expression("GPD log-likelihood for fixed"~beta))

## Plot the profile likelihood for these xi's
## (with an initial interval for the nuisance parameter beta such that
##  logLik_GPD() is finite)
pLL <- sapply(xi., function(xi..) {
    ## Choose beta interval for optimize()
    int <- if(xi.. >= 0) {
               ## Method-of-Moment estimator
               mu.hat <- mean(X)
               sig2.hat <- var(X)
               shape.hat <- (1-mu.hat^2/sig2.hat)/2
               scale.hat <- mu.hat*(1-shape.hat)
               ## log-likelihood always fine for xi.. >= 0 for all beta
               c(1e-8, 2 * scale.hat)
           } else { # xi.. < 0
               ## Make sure logLik_GPD() is finite at endpoints of int
               mx <- max(X)
               -xi.. * mx * c(1.01, 100) # -xi * max(X) * scaling
               ## Note: for shapes xi.. closer to 0, the upper scaling factor
               ##       needs to be chosen sufficiently large in order
               ##       for optimize() to find an optimum (not just the
               ##       upper end point). Try it with '2' instead of '100'.
           }
    ## Optimization
    optimize(function(nuis) logLik_GPD(c(xi.., nuis), x = X),
             interval = int, maximum = TRUE)$maximum
})
plot(xi., pLL, type = "l", xlab = expression(xi),
     ylab = "GPD profile log-likelihood")
}
\keyword{distribution}
