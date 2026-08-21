\name{lcv}
\alias{lcv}
\title{Compute Likelihood Cross Validation Statistic.}
\usage{
lcv(x, \dots)
}
\description{
  The calling sequence for \code{lcv} matches those for the
  \code{\link{locfit}} or \code{\link{locfit.raw}} functions.
  The fit is not returned; instead, the returned object contains
  likelihood cross validation score for the fit.

  The LCV score is exact (up to numerical roundoff) if the
  \code{ev="cross"} argument is provided. Otherwise, the influence
  and cross validated residuals
  are computed using locfit's
  standard interpolation based approximations.
}
\arguments{
  \item{x}{model formula}
  \item{...}{other arguments to locfit}
}

\seealso{
  \code{\link{locfit}},
  \code{\link{locfit.raw}},
  \code{\link{lcvplot}}
}

\keyword{htest}
% Converted by Sd2Rd version 0.2-a5.
