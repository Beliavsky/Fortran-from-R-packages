\name{MultiRNG-package}
\alias{MultiRNG-package}
\alias{MultiRNG}
\docType{package}
\title{
Multivariate Pseudo-Random Number Generation
}
\description{
This package implements the algorithms described in Demirtas (2004) for pseudo-random number generation of 11 multivariate distributions. The following multivariate distributions are available: Normal, \eqn{t}, Uniform, Bernoulli, Hypergeometric, Beta (Dirichlet), Multinomial, Dirichlet-Multinomial, Laplace, Wishart, and Inverted Wishart.

This package contains 11 main functions and 2 auxiliary functions. The methodology for each random-number generation procedure varies and each distribution has its own function. For multivariate normal, \code{draw.d.variate.normal} employs the Cholesky decomposition and a vector of univariate normal draws and for multivariate \eqn{t}, \code{draw.d.variate.t} employs the Cholesky decomposition and a vector of univariate normal and chi-squared draws. \code{draw.d.variate.uniform} is based on cdf of multivariate normal deviates (Falk, 1999) and \code{draw.correlated.binary} generates correlated binary variables using an algorithm developed by Park, Park and Shin (1996) and makes use of the auxiliary function \code{loc.min}. \code{draw.multivariate.hypergeometric} employs sequential generation of succeeding conditionals which are univariate hypergeometric. Furthermore, \code{draw.dirichlet} uses the ratios of gamma variates with a common scale parameter and \code{draw.multinomial} generates data via sequential generation of marginals which are binomials. \code{draw.dirichlet.multinomial} is a mixture distribution of a multinomial that is a realization of a random variable having a Dirichlet distribution. \code{draw.multivariate.laplace} is based on generation of a point s on the d-dimensional sphere and utilizes the auxiliary function \code{generate.point.in.sphere}. \code{draw.wishart} and \code{draw.inv.wishart} employs Wishart variates that follow d-variate normal distribution. 
}
\details{
\tabular{ll}{
Package: \tab MultiRNG\cr
Type: \tab Package\cr
Version: \tab 1.2.4\cr
Date: \tab 2021-03-05\cr
License: \tab GPL-2 | GPL-3\cr
}}
\author{
Hakan Demirtas, Rawan Allozi, Ran Gao

Maintainer: Ran Gao <rgao8@uic.edu>
}
\references{
Demirtas, H. (2004). Pseudo-random number generation in R for commonly used multivariate distributions. \emph{Journal of Modern Applied Statistical Methods}, \bold{3(2)}, 485-497.

Falk, M. (1999). A simple approach to the generation of uniformly distributed random variables with prescribed correlations. \emph{Communications in Statistics, Simulation and Computation}, \bold{28(3)}, 785-791.

Park, C. G., Park, T., & Shin D. W. (1996). A simple method for generating correlated binary variates. \emph{The American Statistician}, \bold{50(4)}, 306-310.
}
