\name{GEV}
\alias{dGEV}
\alias{pGEV}
\alias{qGEV}
\alias{rGEV}
\title{Generalized Extreme Value Distribution}
\description{
  Density, distribution function, quantile function and random variate
  generation for the generalized extreme value distribution (GEV).
}
\usage{
dGEV(x, shape, loc = 0, scale = 1, log = FALSE)
pGEV(q, shape, loc = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
qGEV(p, shape, loc = 0, scale = 1, lower.tail = TRUE, log.p = FALSE)
rGEV(n, shape, loc = 0, scale = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{shape}{GEV shape parameter \eqn{\xi}{xi}, a real.}
  \item{loc}{GEV location parameter \eqn{\mu}{mu}, a real.}
  \item{scale}{GEV scale parameter \eqn{\sigma}{sigma}, a positive real.}
  \item{lower.tail}{\code{\link{logical}}; if \code{TRUE} (default)
    probabilities are \eqn{P(X \le x)}{P(X <= x)} otherwise, \eqn{P(X > x)}.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities \code{p} are
    given as \code{log(p)}.}
}
\value{
  \code{dGEV()} computes the density, \code{pGEV()} the distribution
  function, \code{qGEV()} the quantile function and \code{rGEV()} random
  variates of the generalized extreme value distribution.
}
\details{
  The distribution function of the generalized extreme value
  distribution is given by
  \deqn{F(x) = \left\{ \begin{array}{ll}
    \exp(-(1-\xi(x-\mu)/\sigma)^{-1/\xi}), & \xi\neq 0,\ 1+\xi(x-\mu)/\sigma>0,\\
    \exp(-e^{-(x-\mu)/\sigma}), & \xi = 0,
  \end{array}\right.}{%
    F(x) = exp(-(1-xi(x-mu)/sigma)^{-1/xi}) if xi != 0, 1+xi*(x-mu)/sigma>0 and
    exp(-e^{-(x-mu)/sigma}) if xi = 0,}
  where \eqn{\sigma>0}{sigma>0}.
}
\author{Marius Hofert}
\references{
   McNeil, A. J., Frey, R., and Embrechts, P. (2015).
   \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
   Princeton University Press.
}
\examples{
## Basic sanity checks
plot(pGEV(rGEV(1000, shape = 0.5), shape = 0.5)) # should be U[0,1]
curve(dGEV(x, shape = 0.5), from = -3, to = 5)
}
\keyword{distribution}