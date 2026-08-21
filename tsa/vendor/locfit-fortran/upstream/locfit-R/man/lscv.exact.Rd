\name{lscv.exact}
\alias{lscv.exact}
\title{
Exact LSCV Calculation
}
\usage{
lscv.exact(x, h=0)
}
\description{
  This function performs the exact computation of the least squares
  cross validation statistic for one-dimensional kernel density estimation
  and a constant bandwidth.

  At the time of writing, it is implemented only for the Gaussian
  kernel (with the standard deviation of 0.4; Locfit's standard).
}
\arguments{
  \item{x}{Numeric data vector.}
  \item{h}{The bandwidth. If \code{x} is constructed with
    \code{\link{lp}()}, the bandwidth should be given there instead.}
}
\value{
  A vector of the LSCV statistic and the fitted degrees of freedom.
}
\examples{
data(geyser, package="locfit")
lscv.exact(lp(geyser,h=0.25))
# equivalent form using lscv
lscv(lp(geyser, h=0.25), exact=TRUE)
}
\seealso{
  \code{\link{lscv}},
  \code{\link{lscvplot}}
}

\keyword{htest}
% Converted by Sd2Rd version 0.2-a5.
