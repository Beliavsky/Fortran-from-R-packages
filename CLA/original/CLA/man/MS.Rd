\name{MS}
\alias{MS}
\title{Means (Mu) and Standard Deviations (Sigma) of the \dQuote{Turning Points} from CLA}
\description{
  Compute the vectors of means (\eqn{\mu_i}) and standard deviations
  (\eqn{sigma_i}), for all the turning points of a \code{\link{CLA}}
  result.
}
\usage{
MS(weights_set, mu, covar)
}
\arguments{
  \item{weights_set}{numeric matrix (\eqn{n \times m}{n * m}) of optimal
    asset weights \eqn{W = (w_1, w_2, \ldots, w_m)}, as
    resulting from \code{\link{CLA}()}.}
  \item{mu}{expected (log) returns (identical to \emph{argument} of
    \code{\link{CLA}()}).}
  \item{covar}{covariance matrix of (log) returns (identical to
    \emph{argument} of \code{\link{CLA}()}).}
}
\details{
  These are trivially computable from the \code{\link{CLA}()}'s result.
  To correctly \emph{interpolate} this, \dQuote{hyperbolic}
  interpolation is needed, provided by the \code{\link{findSig}} and
  \code{\link{findMu}} functions.% say more, here or there
}
\value{
  a \code{\link{list}} with components
  \item{Sig}{numeric vector of length \eqn{m} of standard deviations,
    \eqn{\sigma(W)}.}
  \item{Mu}{numeric vector of length \eqn{m} of means \eqn{\mu(W)}.}
}
\author{Yanhao Shi}
\seealso{
 \code{\link{CLA}}.
}
\examples{
## The function is quite simply
MS
## and really an auxiliary function for CLA().

## TODO:  add small (~12 assets) example
}
\keyword{arith}
