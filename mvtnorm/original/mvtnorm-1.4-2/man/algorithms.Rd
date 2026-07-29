\name{algorithms}
\alias{GenzBretz}
\alias{Miwa}
\alias{TVPACK}
\title{ Choice of Algorithm and Hyper Parameters }
\description{
  Choose between three algorithms for evaluating normal (and t-)
  distributions and define hyper parameters.
}
\usage{
GenzBretz(maxpts = 25000, abseps = 0.001, releps = 0)
Miwa(steps = 128, checkCorr = TRUE, maxval = 1e3)
TVPACK(abseps = 1e-6)
}
\arguments{
  \item{maxpts}{maximum number of function values as integer.  The internal
    FORTRAN code always uses a minimum number depending on the dimension.
    (for example 752 for three-dimensional problems).}
  \item{abseps}{absolute error tolerance; for \code{TVPACK} only used
    for dimension 3.}
  \item{releps}{ relative error tolerance as double. }
  \item{steps}{number of grid points to be evaluated; cannot be larger than
    4097.}
  \item{checkCorr}{logical indicating if a check for singularity of the
    correlation matrix should be performed (once per function call to
    \code{pmvt()} or \code{pmvnorm()}).}
  \item{maxval}{replacement for \code{Inf} when non-orthant probabilities
                involving \code{Inf} shall be computed.}
}
\details{

  There are three algorithms available for evaluating normal
  (and two algorithms for t-)
probabilities: The default is the randomized Quasi-Monte-Carlo procedure
by \bibcitet{mvtnorm::numerical-:1992,mvtnorm::comparison:1993} and
\bibcitet{mvtnorm::Genz_Bretz_2002} applicable to
arbitrary covariance structures and dimensions up to 1000.

For normal probabilities, smaller dimensions (up to 20) and non-singular
covariance matrices,
the algorithm by \bibcitet{mvtnorm::Miwa+Hayter+Kuriki:2003} can be used as well. This algorithm can
compute orthant probabilities (\code{lower} being \code{-Inf} or
\code{upper} equal to \code{Inf}). Non-orthant probabilities are computed
from the corresponding orthant probabilities, however, infinite limits are
replaced by \code{maxval} along with a warning.

For two- and three-dimensional problems and semi-infinite integration
region, \code{TVPACK} implements an interface to the methods described
by \bibcitet{mvtnorm::Genz:2004}.

}
\value{
  An object of class \code{"GenzBretz"}, \code{"Miwa"}, or \code{"TVPACK"}
  defining hyper parameters.
}
\references{\bibshow{*}}
\keyword{distribution}
