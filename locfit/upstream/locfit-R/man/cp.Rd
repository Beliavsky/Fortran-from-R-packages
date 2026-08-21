\name{cp}
\alias{cp}
\title{
Compute Mallows' Cp for local regression models.
}
\usage{
cp(x, \dots, sig2=1)
}
\description{
  The calling sequence for \code{cp} matches those for the
  \code{\link{locfit}} or \code{\link{locfit.raw}} functions.
  The fit is not returned; instead, the returned object contains
  Cp criterion for the fit.

  Cp is usually computed using a variance estimate from the largest
  model under consideration, rather than
  \eqn{\sigma^2=1}. This will be done
  automatically when the \code{\link{cpplot}} function is used.

  The Cp score is exact (up to numerical roundoff) if the
  \code{ev="data"} argument is provided. Otherwise, the residual
  sum-of-squares and degrees of freedom are computed using locfit's
  standard interpolation based approximations.
}
\arguments{
  \item{x}{model formula or numeric vector of the independent variable.}
  \item{...}{other arguments to \code{\link{locfit}} and/or
    \code{\link{locfit.raw}}.}
  \item{sig2}{residual variance estimate.}
}

\seealso{
  \code{\link{locfit}},
  \code{\link{locfit.raw}},
  \code{\link{cpplot}}
}

\keyword{htest}
% Converted by Sd2Rd version 0.2-a5.
