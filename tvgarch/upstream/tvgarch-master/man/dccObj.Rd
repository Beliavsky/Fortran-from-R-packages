\name{dccObj}
\alias{dccObj}

\title{Auxiliary functions}

\description{
Auxiliary functions used in the estimation of the multivariate TV(s)-GARCH(p,q,r)-X model. Not intended for the average user.
}

\usage{
  dccObj(par.dcc, z, sigma2, flag)
}

\arguments{
  \item{par.dcc}{numeric \code{vector} containing the ARCH- and GARCH-type coefficients in the dynamic conditional correlations.}
  \item{z}{\code{matrix} of standardized residuals.}
  \item{sigma2}{\code{matrix} of conditional variances.}
  \item{flag}{\code{integer}. If 0, returns a numeric vector with the values of the objective function; if 1 returns the the value of the objective function; if 2, returns the fitted variance components.}
}

\value{
The values of the objective function or fitted dynamic conditional correlations.
}

\author{Susana Campos-Martins}

\references{
Susana Campos-Martins and Genaro Sucarrat (2024) Modeling Nonstationary 
Financial Volatility with the R Package tvgarch, Journal of Statistical 
Software 108, 1-38.

Robert F. Engle (2002) Dynamic conditional correlation: A simple class of 
multivariate generalized autoregressive conditional heteroskedasticity models, 
Journal of Business and Economic Statistics 20, 339-350.
}

\seealso{
\code{\link{mtvgarch}}, 
\code{\link{fitted.mtvgarch}}, 
\code{\link{residuals.mtvgarch}} 
}

\keyword{Time Varying Parameter Models}
\keyword{Nonlinear Time Series}
\keyword{Econometrics}
\keyword{Financial Econometrics}
