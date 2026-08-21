\name{lf}
\alias{lf}
\title{
Locfit term in Additive Model formula
}
\usage{
lf(..., alpha=0.7, deg=2, scale=1, kern="tcub", ev=rbox(), maxk=100)
}
\description{
  This function is used to specify a smooth term in a \code{gam()}
  model formula.

  This function is designed to be used with the S-Plus
  \code{gam()} function. For R users, there are at least two different
  \code{gam()} functions available. Most current distributions of R
  will include the \code{mgcv} library by Simon Wood; \code{lf()}
  is not compatable with this function.

  On CRAN, there is a \code{gam} package by Trevor Hastie, similar to
  the S-Plus version. \code{lf()} should be compatable with this, although
  it's untested.
}
\arguments{
  \item{...}{numeric predictor variable(s)}
  \item{alpha, deg, scale, kern, ev, maxk}{these are as in
    \code{\link{locfit.raw}}.}
}
\examples{
\dontrun{  
# fit an additive semiparametric model to the ethanol data.
stopifnot(require(gam))
# The `gam' package must be attached _before_ `locfit', otherwise
# the following will not work.
data(ethanol, package = "lattice")
fit <- gam(NOx ~ lf(E) + C, data=ethanol)
op <- par(mfrow=c(2, 1))
plot(fit)
par(op)
}
}
\seealso{
  \code{\link{locfit}},
  \code{\link{locfit.raw}},
  \code{\link{gam.lf}},
  \code{gam}
}

\keyword{models}
