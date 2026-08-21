\name{lscv}
\alias{lscv}
\title{
Least Squares Cross Validation Statistic.
}
\usage{
lscv(x, \dots, exact=FALSE)
}
\description{
  The calling sequence for \code{lscv} matches those for the
  \code{\link{locfit}} or \code{\link{locfit.raw}} functions.
  Note that this function is only designed for density estimation
  in one dimension. The returned object contains the
  least squares cross validation score for the fit.

  The computation of \eqn{\int \hat f(x)^2 dx} is performed numerically.
  For kernel density estimation, this is unlikely to agree exactly
  with other LSCV routines, which may perform the integration analytically.
}
\arguments{
  \item{x}{model formula (or numeric vector, if \code{exact=T})}
  \item{...}{other arguments to \code{\link{locfit}} or
    \code{\link{lscv.exact}} }
  \item{exact}{By default, the computation is approximate.
    If \code{exact=TRUE}, exact computation using
    \code{\link{lscv.exact}} is performed. This uses kernel density estimation
    with a constant bandwidth.}
}
\value{
  A vector consisting of the LSCV statistic and fitted degrees of freedom.
}
\examples{
# approximate calculation for a kernel density estimate
data(geyser, package="locfit")
lscv(~lp(geyser,h=1,deg=0), ev=lfgrid(100,ll=1,ur=6), kern="gauss")
# same computation, exact
lscv(lp(geyser,h=1),exact=TRUE)
}
\seealso{
  \code{\link{locfit}},
  \code{\link{locfit.raw}},
  \code{\link{lscv.exact}}
  \code{\link{lscvplot}}
}

\keyword{htest}
% Converted by Sd2Rd version 0.2-a5.
