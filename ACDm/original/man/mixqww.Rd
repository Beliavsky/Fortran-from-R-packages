\name{Discreetly mixed Q-Weibull and ordinary Weibull}
\alias{Discreetly mixed Q-Weibull and ordinary Weibull}
\alias{dmixqww}
\alias{pmixqww}
\alias{mixqwwHazard}

\title{
Discreet mix of the q-Weibull and the ordinary Weibull distributions
}
\description{
Density (PDF), distribution function (CDF), and hazard function for a discreetly mixed distribution of the Q-Weibull and the ordinary Weibull distributions.
}
\usage{
dmixqww(x, pdist = .5, a = 1.2, qdist = 1.5, theta = .8, gamma = 1, b = 1,
  forceExpectation = F)
        
pmixqww(q, pdist = .5, a = 1.2, qdist = 1.5, theta = .8, gamma = 1, b = 1,
  forceExpectation = F)
        
mixqwwHazard(x, pdist = .5, a = 1.2, qdist = 1.5, theta = .8, gamma = 1, b = 1,
  forceExpectation = F)
}
\arguments{
  \item{x, q}{
vector of quantiles.
}
  \item{pdist, a, qdist, theta, gamma, b}{
parameters, see 'Details'.
}
  \item{forceExpectation}{
logical; if \code{TRUE}, the expectation of the distribution is forced to be 1 by letting \code{b} be a function of the other parameters.
}
}
\details{
The PDF for the mixed distribution is:

\deqn{f(x) = p(2-q)\frac{a}{b^a} x^{a-1} \left[1-(1-q)\left(\frac{x}{b}\right)^a\right]^{\frac{1}{1-q}} + (1-p)\theta \gamma x^{-\theta x^{\gamma}}}

if \code{forceExpectation = TRUE} the \code{b} parameter is a function of the other parameters to force the expectation to be 1.

}

\seealso{
  \code{\link{qWeibullDist}} for the Q-Weibull distribution and \code{\link{pmixqwe}} for Q-Weibull mixed with the exponential distribution.
}