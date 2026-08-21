\name{regband}
\alias{regband}
\title{
Bandwidth selectors for local regression.
}
\usage{
regband(formula, what = c("CP", "GCV", "GKK", "RSW"), deg=1, ...)
}
\description{
Function to compute local regression bandwidths for local linear regression,
implemented as a front end to \code{\link{locfit}()}.

  This function is included for comparative purposes only. Plug-in
selectors are based on flawed logic, make unreasonable and restrictive
assumptions
and do not use the full power of the estimates available in Locfit.
Any relation between the results produced by this function and
desirable estimates are entirely coincidental.
}
\arguments{
\item{formula}{Model Formula (one predictor).}
\item{what}{Methods to use.}
\item{deg}{Degree of fit.}
\item{...}{Other Locfit options.}
}
\value{
Vector of selected bandwidths.
}
\keyword{htest}
