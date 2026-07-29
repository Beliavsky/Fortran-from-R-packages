\name{GeneralizedGammaDist}
\alias{GeneralizedGammaDist}
\alias{dgengamma}
\alias{pgengamma}
\alias{qgengamma}
\alias{rgengamma}
\alias{gengammaHazard}

\title{
The generelized Gamma distribution
}
\description{
Density (PDF), distribution function (CDF), quantile function (inverted CDF), random generation and hazard function for the generelized Gamma distribution with parameters gamma, kappa and lambda.
}
\usage{
dgengamma(x, gamma = 0.3, kappa = 1.2, lambda = 0.3, forceExpectation = F)
pgengamma(x, gamma = .3, kappa = 3, lambda = .3, forceExpectation = F)
qgengamma(p, gamma = .3, kappa = 3, lambda = .3, forceExpectation = F)
rgengamma(n = 1, gamma = .3, kappa = 3, lambda = .3, forceExpectation = F)
gengammaHazard(x, gamma = .3, kappa = 3, lambda = .3, forceExpectation = F)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
vector of quantiles.
}
  \item{p}{
vector of probabilities.
}
  \item{n}{
number of observations..
}
  \item{gamma, kappa, lambda}{
parameters, see 'Details'.
}
  \item{forceExpectation}{
logical; if \code{TRUE}, the expectation of the distribution is forced to be 1 by letting \code{theta} be a function of the other parameters. 
}
}
\details{
The PDF for the generelized Gamma distribution is:

\deqn{f(x)=\frac{\gamma x^{\kappa \gamma - 1}}{\lambda^{\kappa \gamma}\Gamma (\kappa)}\exp \left\{{-\left(\frac{x}{\lambda}\right)^{\gamma}}\right\}}


}
\value{
\code{dgengamma} gives the density (PDF), \code{pgengamma} gives the distribution function (CDF), \code{qgengamma} gives the quantile function (inverted CDF), \code{rgenGamma} generates random deviates, and \code{genGammaHazard} gives the hazard function.

}
\author{
Markus Belfrage
}