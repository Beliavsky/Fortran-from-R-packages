\name{qc}
\alias{qc}
\alias{qc}
\docType{package}
\title{Quantile Correlation
}
\description{
qc is used to compute  the quantile correlation with given quantiles.
}
\usage{
qc(x, y, tau)
}

\arguments{
\item{x}{The covariate variable.}
\item{y}{The response variable.}
\item{tau}{The quantile(s) to be estimated.}
}

\value{
\item{tau}{The quantile(s).}
\item{rho}{The value of quantile correlation.}
}

\author{
\packageAuthor{QCSIS}
}
\references{
 Li et al.. Quantile correlations and quantile autoregressive modeling. Journal of the American Statistical Association,2015,110(509):246--261.\cr
}

\examples{
n   <- 1000
x   <- rnorm(n)
y   <- 2 * x + rt(n,df = 1)
tau <- 1:9 / 10
qc(x = x, y = y, tau = tau)
}
