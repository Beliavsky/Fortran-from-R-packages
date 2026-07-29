\name{riskmeasures}
\alias{riskmeasures}
\title{Calculates Quantiles and Expected Shortfalls}
\description{
  Makes a rapid calculation of point estimates of prescribed
  quantiles and expected shortfalls using the output of the
  function \code{gpd}. 
}
\usage{
riskmeasures(x, p)
}
\arguments{
\item{x}{results of a \code{gpd} fit}
\item{p}{a vector of probability levels}
}
\value{
  A matrix with three columns: probability level, quantile
  estimate, shortfall estimate. 
}
\details{
  This function simply calculates point estimates and (at present)
  makes no attempt to calculate confidence intervals for the risk
  measures. If confidence levels are required use \code{gpd.q} and
  \code{gpd.sfall} which interact with graphs of the tail of a loss
  distribution and are much slower. 
}
\seealso{\code{\link{gpd}}, \code{\link{tailplot}},
  \code{\link{gpd.q}}, \code{\link{gpd.sfall}}}

\examples{
data(danish)
out <- gpd(danish, 10) 
riskmeasures(out, c(0.999, 0.9999)) 
# gives estimates of 0.999 and 0.9999 quantiles of Danish loss
# distribution as well as the associated expected shortfall estimates 
}
\keyword{htest}
