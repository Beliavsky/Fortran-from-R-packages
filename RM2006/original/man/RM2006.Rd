\name{RM2006}
\alias{RM2006}
\title{
RiskMetrics 2006 Methodology
}
\description{
Estimation of the conditional covariance matrix using the RiskMetrics 2006 methodology of Zumbach (2007).
}
\usage{
RM2006(data, tau0, tau1, kmax,rho)
}
\arguments{
  \item{data}{
Matrix containing a TxK time series returns.
}
  \item{tau0}{
optional input parameter. Default 1560
}
  \item{tau1}{
optional input parameter. Default 4
}
  \item{kmax}{
optional input parameter. Default 14
}
  \item{rho}{
optional input parameter. Default 1.4142
}
}
\details{
More details can be found in Zumbach (2007) and in the MFE Toolbox of Kevin Sheppard (function riskmetrics2006).
}
\value{
The funcion returns an array containing for each t (t = 1, ..., T+1) a KxK matrix with the conditional covariance matrix estimates.
}
\references{
Zumbach, G. (2007) The Riskmetrics 2006 Methodology. Available at SSRN: https://ssrn.com/abstract=1420185 or http://dx.doi.org/10.2139/ssrn.1420185
}
\author{
Carlos Trucios
}
\examples{
Data=matrix(rnorm(1000),nrow = 100, ncol = 10)
RM2006(Data)
}

\keyword{Risk Metrics}
\keyword{Conditional Covariance Matrix}
