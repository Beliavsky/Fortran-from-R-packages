\name{SkTDist}
\alias{SkTDist}
\alias{dskt}
\alias{pskt}
\alias{qskt}
\alias{rskt}
\title{The Skewed Student t Distribution}
\description{
  Density, distribution function, quantile function and random
  generation for the skewed t distribution, as introduced by Fernandez and
Steel, with \code{df} degrees of freedom.
}
\usage{
dskt(x, df, gamma = 1)
pskt(x, df, gamma = 1)
qskt(p, df, gamma)
rskt(n, df, gamma)
}
\arguments{
  \item{x}{vector of quantiles.}
  \item{p}{vector of probabilities.}
  \item{n}{number of observations. If \code{length(n) > 1}, the length
    is taken to be the number required.}
  \item{df}{degrees of freedom (\eqn{> 0}, maybe non-integer).}
  \item{gamma}{skewing parameter, \eqn{\gamma}{gamma}}
}
\value{
  \code{dskt} gives the density,
  \code{pskt} gives the distribution function,
  \code{qskt} gives the quantile function, and
  \code{rskt} generates random deviates.
}
\details{
  The Skewed \eqn{t} distribution with \code{df} \eqn{= \nu}{= n} degrees of
  freedom has the following density, where \eqn{f(x)} is the density of the
\eqn{t} distribution, with  \eqn{= \nu}{= n} degrees of
  freedom :
\deqn{f(x) = \frac{2}{\gamma + \frac{1}{\gamma}} f(\gamma x) \quad for 
\quad x<0}{f(x) = 2/(gamma + 1/gamma) f(x gamma) for x < 0} and
\deqn{f(x) = \frac{2}{\gamma + \frac{1}{\gamma}} f(\frac{x}{\gamma}) \quad
  for \quad x \ge 0}{f(x) = 2/(gamma + 1/gamma) f(x/gamma) for x <= 0}

}
\references{
  Fernandez, C. and Steel, M. F. J. (1998).
  On Bayesian modeling of fat tails and skewness,
  \emph{J. Am. Statist. Assoc.} \bold{93}, 359--371.

  Rohr, P. and Hoeschele, I. (2002).
  Bayesian QTL mapping using skewed Student-\eqn{t} distributions,
  \emph{Genet. Sel. Evol.} \bold{34}, 1--21.
}
\seealso{\code{\link{df}} for the F distribution.}
\examples{
dskt(0.5,2)
dskt(0.01,2,2)
pskt(1.25,2,2)
pskt(c(0.5,1.25),3)
qskt(c(0,0.025,0.25,0.5,0.75,0.975,1),2,2)
rskt(100,2,2)
plot(function(x)dskt(x,2,2),-3,3,n=301)
}
\keyword{distribution}
