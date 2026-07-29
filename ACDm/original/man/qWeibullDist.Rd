\name{qWeibullDist}
\alias{qWeibullDist}
\alias{dqweibull}
\alias{pqweibull}
\alias{qqweibull}
\alias{rqweibull}
\alias{qweibullExpectation}
\alias{qweibullHazard}

\title{
The q-Weibull distribution
}
\description{
Density (PDF), distribution function (CDF), quantile function (inverted CDF), random generation, exepcted value, and hazard function for the q-Weibull distribution.
}
\usage{
dqweibull(x, a = .8, qdist = 1.2, b = 1, forceExpectation = F)
pqweibull(q, a = .8, qdist = 1.2, b = 1, forceExpectation = F)
qqweibull(p, a = .8, qdist = 1.2, b = 1, forceExpectation = F)
rqweibull(n = 1, a = .8, qdist = 1.2, b = 1, forceExpectation = F)
qweibullExpectation(a = .8, qdist = 1.2, b = 1)
qweibullHazard(x, a = .8, qdist = 1.2, b = 1, forceExpectation = F)
}
\arguments{
  \item{x, q}{
vector of quantiles.
}
  \item{p}{
vector of probabilities.
}
  \item{n}{
number of observations.
}
  \item{a, qdist, b}{
parameters, see 'Details'.
}
  \item{forceExpectation}{
logical; if \code{TRUE}, the expectation of the distribution is forced to be 1 by letting \code{b} be a function of the other parameters.
}
}
\details{
The PDF for the q-Weibull distribution is:

\deqn{f(\epsilon) = (2-q)\frac{a}{b^a} \epsilon^{a-1} \left[1-(1-q)\left(\frac{\epsilon}{b}\right)^a\right]^{\frac{1}{1-q}}}

The distribution was used for ACD models by Vuorenmaa (2009).
}

\references{
Vuorenmaa, T. (2009) 
\emph{A q-Weibull Autoregressive Conditional Duration Model with an Application to NYSE and HSE data.} Available at SSRN: http://ssrn.com/abstract=1952550.

}