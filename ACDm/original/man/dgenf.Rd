\name{dgenf}
\alias{dgenf}
\alias{pgenf}
\alias{genfHazard}

\title{
The generalized F distribution
}
\description{
Density and distribution function for the generalized F distribution. Warning: the distribution function \code{pgenf} and \code{genfHazard} are computed numerically, and may not be precise!
}
\usage{
dgenf(x, kappa = 5, eta = 1.5, gamma = .8, lambda = 1, forceExpectation = F)
pgenf(q, kappa = 5, eta = 1.5, gamma = .8, lambda = 1, forceExpectation = F)
genfHazard(x, kappa = 5, eta = 1.5, gamma = .8, lambda = 1, forceExpectation = F)
}
\arguments{
  \item{x, q}{
vector of quantiles.
}
  \item{kappa, eta, gamma, lambda}{
parameters, see 'Details'.
}
  \item{forceExpectation}{
logical; if \code{TRUE}, the expectation of the distribution is forced to be 1 by letting theta be a function of the other parameters.
}
}
\details{
The PDF for the generelized F distribution is:
\deqn{f(\epsilon)= \frac{\gamma \epsilon^{\kappa \gamma -1}[\eta+(\epsilon/\lambda)^{\gamma}]^{-\eta-\kappa}\eta^{\eta}}{\lambda^{\kappa \gamma}B(\kappa,\eta)},}
where \eqn{B(\kappa,\eta)=\frac{\Gamma(\kappa)\Gamma(\eta)}{\Gamma(\kappa+\eta)}} is the beta function.
}