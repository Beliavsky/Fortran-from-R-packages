\name{kdeb}
\alias{kdeb}
\title{
Bandwidth selectors for kernel density estimation.
}
\usage{
kdeb(x, h0 = 0.01 * sd, h1 = sd, meth = c("AIC", "LCV", "LSCV", "BCV", 
  "SJPI", "GKK"), kern = "gauss", gf = 2.5)

}
\description{
Function to compute kernel density estimate bandwidths, as used in the
simulation results in Chapter 10 of Loader (1999).

  This function is included for comparative purposes only. Plug-in
selectors are based on flawed logic, make unreasonable and restrictive
assumptions
and do not use the full power of the estimates available in Locfit.
Any relation between the results produced by this function and
desirable estimates are entirely coincidental.
}
\arguments{
\item{x}{One dimensional data vector.}
\item{h0}{Lower limit for bandwidth selection. Can be fairly small,
  but h0=0 would cause problems.}
\item{h1}{Upper limit.}
\item{meth}{Required selection method(s).}
\item{kern}{Kernel. Most methods require \code{kern="gauss"}, the default
  for this function only.}
\item{gf}{Standard deviation for the gaussian kernel. Default 2.5, as
  Locfit's standard. Most papers use 1.
}
}
\value{
Vector of selected bandwidths.
}
\references{
Loader, C. (1999). Local Regression and Likelihood. Springer, New York.
}
\keyword{htest}
