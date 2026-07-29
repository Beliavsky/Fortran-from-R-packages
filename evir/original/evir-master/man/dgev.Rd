\name{dgev}
\alias{dgev}
\alias{pgev}
\alias{qgev}
\alias{rgev}
\title{Generalized Extreme Value Distribution}
\description{
  Cumulative probability, quantiles, density and random generation
  from the generalized extreme value distribution. 
}
\usage{
pgev(q, xi = 1, mu = 0, sigma = 1) 
qgev(p, xi = 1, mu = 0, sigma = 1) 
dgev(x, xi = 1, mu = 0, sigma = 1) 
rgev(n, xi = 1, mu = 0, sigma = 1)
}
\arguments{
\item{q}{vector of quantiles}
\item{p}{vector of probabilities}
\item{x}{vector of values at which to evaluate density}
\item{n}{sample size}
\item{xi}{shape parameter}
\item{mu}{location parameter}
\item{sigma}{scale parameter}
}
\value{
  Probability (\code{pgev}), quantile (\code{qgev}), density
  (\code{dgev}) or random sample (\code{rgev}) for the GEV
  distribution with shape \code{xi}. 
}

\seealso{\code{\link{dgpd}}, \code{\link{gev}}}
\keyword{distribution}
