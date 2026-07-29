\name{BurrDist}
\alias{BurrDist}
\alias{dburr}
\alias{pburr}
\alias{qburr}
\alias{rburr}
\alias{burrExpectation}

\title{
The Burr Distribution
}

\description{
Density, distribution function, quantile function, random generation and calculation of the expected value for the Burr distribution with parameters theta, kappa and sig2.
}

\usage{
dburr(x, theta = 1, kappa = 1.2, sig2 = 0.3, forceExpectation = F)
pburr(x, theta = 1, kappa = 1.2, sig2 = .3, forceExpectation = F)
qburr(p, theta = 1, kappa = 1.2, sig2 = .3, forceExpectation = F)
rburr(n = 1, theta = 1, kappa = 1.2, sig2 = .3, forceExpectation = F)
burrExpectation(theta = 1, kappa = 1.2, sig2 = .3)
}


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
  \item{theta, kappa, sig2}{
parameters, see 'Details'.
}
  \item{forceExpectation}{
logical; if \code{TRUE}, the expectation of the distribution is forced to be 1 by letting \code{theta} be a function of the other parameters. 
}
}
\details{
The PDF for the Burr distribution is (as in e.g. Grammig and Maurer, 2000):
\deqn{f(x)=\frac{\theta \kappa x^{\kappa - 1}}{(1 + \sigma^2 x^{\kappa)^{\frac{1}{\sigma^2}+1}}}}

}
\value{

\code{dburr} gives the density (PDF), \code{qburr} the quantile function (inverted CDF), \code{rburr} generates random deviates, and \code{burrExpectation} returns the expected value of the distribution, given the parameters.

}
\references{

Grammig, J., and Maurer, K.-O. (2000)
\emph{Non-monotonic hazard functions and the autoregressive conditional duration model}.
Econometrics Journal 3: 16-38.

}
\author{
Markus Belfrage
}