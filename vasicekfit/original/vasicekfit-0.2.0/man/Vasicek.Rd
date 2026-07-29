\name{Vasicek}
\alias{dvasicek}
\alias{pvasicek}
\alias{qvasicek}
\alias{rvasicek}
\title{The Vasicek Loss Distribution}
\description{
Density, distribution function, quantile function, and random generation
for the Vasicek credit loss distribution.
}
\usage{
dvasicek(x, p, rho, kappa = NULL, u = NULL, log = FALSE)
pvasicek(q, p, rho, kappa = NULL, u = NULL, lower.tail = TRUE, log.p = FALSE)
qvasicek(prob, p, rho, kappa = NULL, u = NULL, lower.tail = TRUE, log.p = FALSE)
rvasicek(n, p, rho, kappa = NULL, u = NULL)
}
\arguments{
  \item{x, q}{vector of quantiles (loss rates in (0, 1)).}
  \item{prob}{vector of probabilities.}
  \item{n}{number of observations to generate.}
  \item{p}{probability of default, in (0, 1).}
  \item{rho}{asset value correlation, in (0, 1).}
  \item{kappa}{optional numeric vector of macro-factor sensitivities.}
  \item{u}{optional numeric vector of macro-factor values, same length as
    \code{kappa}.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities/densities are given
    as \code{log(p)}.}
  \item{lower.tail}{logical; if \code{TRUE} (default), probabilities are
    \eqn{P[X \le x]}{P[X <= x]}, otherwise \eqn{P[X > x]}.}
}
\details{
The Vasicek loss distribution arises from the single-factor Gaussian copula
model. Its probability density function is

\deqn{f(x) = \sqrt{\frac{1-\rho}{\rho}} \exp\left\{\frac{1}{2}\left(\Phi^{-1}(x)^2 - \left[\frac{\sqrt{1-\rho}\,\Phi^{-1}(x) - \sum \kappa_j u_j - \Phi^{-1}(p)}{\sqrt{\rho}}\right]^2\right)\right\}}{f(x) = sqrt((1-rho)/rho) * exp(0.5 * (qnorm(x)^2 - ((sqrt(1-rho)*qnorm(x) - sum(kappa * u) - qnorm(p))/sqrt(rho))^2))}

When \code{kappa} and \code{u} are \code{NULL}, this reduces to the standard
two-parameter Vasicek distribution.
}
\value{
\code{dvasicek} gives the density, \code{pvasicek} gives the distribution
function, \code{qvasicek} gives the quantile function, and \code{rvasicek}
generates random deviates.
}
\references{
Vasicek, O. A. (2002). The distribution of loan portfolio value.
\emph{Risk}, 15(12), 160--162.
}
\examples{
# Standard Vasicek density
curve(dvasicek(x, p = 0.03, rho = 0.10), from = 0.001, to = 0.20,
      ylab = "Density", main = "Vasicek loss distribution")

# 99th percentile (Basel-style)
qvasicek(0.99, p = 0.03, rho = 0.10)

# Check CDF and quantile are inverses
qvasicek(pvasicek(0.05, p = 0.03, rho = 0.10), p = 0.03, rho = 0.10)

# Random sample
set.seed(1)
hist(rvasicek(10000, p = 0.03, rho = 0.10), breaks = 100,
     main = "Vasicek loss distribution samples")
}
\seealso{
\code{\link{vasicekfit}}
}
