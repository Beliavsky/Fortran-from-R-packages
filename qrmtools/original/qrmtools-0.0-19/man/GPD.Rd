\name{GPD}
\alias{dGPD}
\alias{pGPD}
\alias{qGPD}
\alias{rGPD}
\alias{dPar}
\alias{pPar}
\alias{qPar}
\alias{rPar}
\title{(Generalized) Pareto Distribution}
\description{
  Density, distribution function, quantile function and random variate
  generation for the (generalized) Pareto distribution (GPD).
}
\usage{
dGPD(x, shape, scale, log = FALSE)
pGPD(q, shape, scale, lower.tail = TRUE, log.p = FALSE)
qGPD(p, shape, scale, lower.tail = TRUE, log.p = FALSE)
rGPD(n, shape, scale)

dPar(x, shape, scale = 1, log = FALSE)
pPar(q, shape, scale = 1, lower.tail = TRUE, log.p = FALSE)
qPar(p, shape, scale = 1, lower.tail = TRUE, log.p = FALSE)
rPar(n, shape, scale = 1)
}
\arguments{
  \item{x, q}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations.}
  \item{shape}{GPD shape parameter \eqn{\xi}{xi} (a real number) and
    Pareto shape parameter \eqn{\theta}{theta} (a positive number).}
  \item{scale}{GPD scale parameter \eqn{\beta}{beta} (a positive number)
    and Pareto scale parameter \eqn{\kappa}{kappa} (a positive number).}
  \item{lower.tail}{\code{\link{logical}}; if \code{TRUE} (default)
    probabilities are \eqn{P(X \le x)}{P(X <= x)} otherwise, \eqn{P(X > x)}.}
  \item{log, log.p}{logical; if \code{TRUE}, probabilities \code{p} are
    given as \code{log(p)}.}
}
\value{
  \code{dGPD()} computes the density, \code{pGPD()} the distribution
  function, \code{qGPD()} the quantile function and \code{rGPD()} random
  variates of the generalized Pareto distribution.

  Similary for \code{dPar()}, \code{pPar()}, \code{qPar()} and
  \code{rPar()} for the Pareto distribution.
}
\details{
  The distribution function of the generalized Pareto distribution is
  given by
  \deqn{F(x) = \left\{ \begin{array}{ll}
    1-(1+\xi x/\beta)^{-1/\xi}, & \xi \neq 0,\\
    1-\exp(-x/\beta), & \xi = 0,
    \end{array}\right.}{%
    F(x) = 1-(1+xi x/beta)^{-1/xi} if xi != 0 and
    1-exp(-x/beta) if xi = 0,}
  where \eqn{\beta>0}{beta>0} and \eqn{x\ge0}{x >= 0} if \eqn{\xi\ge
  0}{xi >= 0}
  and \eqn{x\in[0,-\beta/\xi]}{x in [0,-beta/xi]} if \eqn{\xi<0}{xi<0}.

  The distribution function of the Pareto distribution is given by
  \deqn{F(x) = 1-(1+x/\kappa)^{-\theta},\ x\ge 0,}{F(x) = 1-(1+x/kappa)^{-theta}, x
  >= 0,} where \eqn{\theta > 0}{theta > 0}, \eqn{\kappa > 0}{kappa > 0}.

  In contrast to \code{dGPD()}, \code{pGPD()}, \code{qGPD()} and
  \code{rGPD()}, the functions \code{dPar()}, \code{pPar()},
  \code{qPar()} and \code{rPar()} are vectorized in their main
  argument and the parameters.
}
\author{Marius Hofert}
\references{
   McNeil, A. J., Frey, R., and Embrechts, P. (2015).
   \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
   Princeton University Press.
}
\examples{
## Basic sanity checks
curve(dGPD(x, shape = 0.5, scale = 3), from = -1, to = 5)
plot(pGPD(rGPD(1000, shape = 0.5, scale = 3), shape = 0.5, scale = 3)) # should be U[0,1]
}
\keyword{distribution}