\name{Discreetly mixed Q-Weibull and exponential}
\alias{Discreetly mixed Q-Weibull and exponential}
\alias{dmixqwe}
\alias{pmixqwe}
\alias{mixqweHazard}

\title{
Discreet mix of the Q-Weibull and the exponential distributions
}
\description{
Density (PDF), distribution function (CDF), and hazard function for a discreetly mixed distribution of the Q-Weibull and the exponential distributions.
}
\usage{
dmixqwe(x, pdist = .5, a = .8, qdist = 1.5, lambda = .8, b = 1, forceExpectation = F)
pmixqwe(q, pdist = .5, a = .8, qdist = 1.5, lambda = .8, b = 1, forceExpectation = F)
mixqweHazard(x, pdist = .5, a = .8, qdist = 1.5, lambda = .8, b = 1, forceExpectation = F)
}
\arguments{
  \item{x, q}{
vector of quantiles.
}
  \item{pdist, a, qdist, lambda, b}{
parameters, see 'Details'.
}
  \item{forceExpectation}{
logical; if \code{TRUE}, the expectation of the distribution is forced to be 1 by letting \code{b} be a function of the other parameters.
}
}
\details{
The PDF for the mixed distribution is:

\deqn{f(x) = p(2-q)\frac{a}{b^a} x^{a-1} \left[1-(1-q)\left(\frac{x}{b}\right)^a\right]^{\frac{1}{1-q}} + (1-p)\frac{1}{\lambda}exp(-\frac{x}{\lambda})}

if \code{forceExpectation = TRUE} the \code{b} parameter is a function of the other parameters to force the expectation to be 1.

}

\seealso{
  \code{\link{qWeibullDist}} for the Q-Weibull distribution and \code{\link{pmixqww}} for Q-Weibull mixed with the ordinary Weibull.
}