\name{dgpd}
\alias{dgpd}
\alias{pgpd}
\alias{qgpd}
\alias{rgpd}
\title{Generalized Pareto Distribution }
\description{
  Cumulative probability, quantiles, density and random
  generation from the generalized Pareto distribution. 
}
\usage{
pgpd(q, xi, mu = 0, beta = 1) 
qgpd(p, xi, mu = 0, beta = 1) 
dgpd(x, xi, mu = 0, beta = 1) 
rgpd(n, xi, mu = 0, beta = 1)
}
\arguments{
\item{q}{vector of quantiles.}
\item{p}{vector of probabilities.}
\item{x}{vector of values at which to evaluate density}
\item{n}{sample size}
\item{xi}{shape parameter.}
\item{mu}{location parameter.}
\item{beta}{scale parameter}
}
\value{
  Probability (\code{pgpd}), quantile (\code{qgpd}), density
  (\code{dgpd}) or random sample (\code{rgpd}) for the GPD
  distribution with shape \code{xi}. 
}

\seealso{\code{\link{dgev}}, \code{\link{gpd}}}
\keyword{distribution}
