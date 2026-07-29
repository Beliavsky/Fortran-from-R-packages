\name{GPDtail}
\alias{dGPDtail}
\alias{pGPDtail}
\alias{qGPDtail}
\alias{rGPDtail}
\title{GPD-Based Tail Distribution (POT method)}
\description{
  Density, distribution function, quantile function and random variate
  generation for the GPD-based tail distribution in the POT method.
}
\usage{
dGPDtail(x, threshold, p.exceed, shape, scale, log = FALSE)
pGPDtail(q, threshold, p.exceed, shape, scale, lower.tail = TRUE, log.p = FALSE)
qGPDtail(p, threshold, p.exceed, shape, scale, lower.tail = TRUE, log.p = FALSE)
rGPDtail(n, threshold, p.exceed, shape, scale)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{threshold}{threshold \eqn{u} in the POT method.}
  \item{p.exceed}{probability of exceeding the threshold u; for the
    Smith estimator, this is \code{mean(x > threshold)} for \code{x}
    being the data.}
  \item{shape}{GPD shape parameter \eqn{\xi}{xi} (a real number).}
  \item{scale}{GPD scale parameter \eqn{\beta}{beta} (a positive number).}
  \item{lower.tail}{\code{\link{logical}}; if \code{TRUE} (default)
    probabilities are \eqn{P(X \le x)}{P(X <= x)} otherwise, \eqn{P(X > x)}.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities \code{p} are
    given as \code{log(p)}.}
}
\value{
  \code{dGPDtail()} computes the density, \code{pGPDtail()} the distribution
  function, \code{qGPDtail()} the quantile function and \code{rGPDtail()} random
  variates of the GPD-based tail distribution in the POT method.
}
\details{
  Let \eqn{u} denote the threshold (\code{threshold}), \eqn{p_u} the exceedance
  probability (\code{p.exceed}) and \eqn{F_{GPD}} the GPD
  distribution function. Then the distribution function of the GPD-based tail
  distribution is given by
  \deqn{F(q) = 1-p_u(1-F_{GPD}(q - u))}. The quantile function is
  \deqn{F^{-1}(p) = u + F_GPD^{-1}(1-(1-p)/p_u)} and the density is
  \deqn{f(x) = p_u f_{GPD}(x - u)}, where \eqn{f_{GPD}} denotes the GPD
  density.

  Note that the distribution function has a jumpt of height \eqn{P(X \le
    u)}{P(X <=u)} (\code{1-p.exceed}) at \eqn{u}.
}
\author{Marius Hofert}
\references{
   McNeil, A. J., Frey, R., and Embrechts, P. (2015).
   \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
   Princeton University Press.
}
\examples{
## Generate data to work with
set.seed(271)
X <- rt(1000, df = 3.5) # in MDA(H_{1/df}); see MFE (2015, Section 16.1.1)

## Determine thresholds for POT method
mean_excess_plot(X[X > 0])
abline(v = 1.5)
u <- 1.5 # threshold

## Fit GPD to the excesses (per margin)
fit <- fit_GPD_MLE(X[X > u] - u)
fit$par
1/fit$par["shape"] # => close to df

## Estimate threshold exceedance probabilities
p.exceed <- mean(X > u)

## Define corresponding densities, distribution function and RNG
dF <- function(x) dGPDtail(x, threshold = u, p.exceed = p.exceed,
                           shape = fit$par["shape"], scale = fit$par["scale"])
pF <- function(q) pGPDtail(q, threshold = u, p.exceed = p.exceed,
                           shape = fit$par["shape"], scale = fit$par["scale"])
rF <- function(n) rGPDtail(n, threshold = u, p.exceed = p.exceed,
                           shape = fit$par["shape"], scale = fit$par["scale"])

## Basic check of dF()
curve(dF, from = u - 1, to = u + 5)

## Basic check of pF()
curve(pF, from = u, to = u + 5, ylim = 0:1) # quite flat here
abline(v = u, h = 1-p.exceed, lty = 2) # mass at u is 1-p.exceed (see 'Details')

## Basic check of rF()
set.seed(271)
X. <- rF(1000)
plot(X., ylab = "Losses generated from the fitted GPD-based tail distribution")
stopifnot(all.equal(mean(X. == u), 1-p.exceed, tol = 7e-3)) # confirms the above
## Pick out 'continuous part'
X.. <- X.[X. > u]
plot(pF(X..), ylab = "Probability-transformed tail losses") # should be U[1-p.exceed, 1]
}
\keyword{distribution}