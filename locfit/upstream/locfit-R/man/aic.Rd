\name{aic}
\alias{aic}
\title{
Compute Akaike's Information Criterion.
}
\usage{
aic(x, \dots, pen=2)
}
\description{
  The calling sequence for \code{aic} matches those for the
  \code{\link{locfit}} or \code{\link{locfit.raw}} functions.
  The fit is not returned; instead, the returned object contains
  Akaike's information criterion for the fit.

  The definition of AIC used here is -2*log-likelihood + pen*(fitted d.f.).
  For quasi-likelihood, and local regression, this assumes the scale
  parameter is one. Other scale parameters can effectively be used
  by changing the penalty.

  The AIC score is exact (up to numerical roundoff) if the
  \code{ev="data"} argument is provided. Otherwise, the residual
  sum-of-squares and degrees of freedom are computed using locfit's
  standard interpolation based approximations.
}
\arguments{
  \item{x}{model formula}
  \item{...}{other arguments to locfit}
  \item{pen}{penalty for the degrees of freedom term}
}

\seealso{
  \code{\link{locfit}},
  \code{\link{locfit.raw}},
  \code{\link{aicplot}}
}

\keyword{htest}
% Converted by Sd2Rd version 0.2-a5.
