\name{maxLik-package}
\alias{maxLik-package}
\docType{package}
\title{
Maximum Likelihood Estimation
}
\description{
  This package contains a set of functions and tools for Maximum Likelihood (ML)
  estimation.  The focus of the package is on non-linear
  optimization from the ML viewpoint, and it provides several convenience wrappers
  and tools, like BHHH algorithm, variance-covariance
  matrix and standard errors.
}
\details{
\pkg{maxLik} package is a set of convenience tools and wrappers
focusing on
Maximum Likelihood (ML) analysis, but it also contains tools for
other optimization tasks.
The package includes a) wrappers for several
existing optimizers (implemented by \code{\link[stats:optim]{optim}}); b) original
optimizers, including Newton-Raphson and Stochastic Gradient Ascent;
and c) several convenience tools
to use these optimizers from the ML perspective.  Examples are BHHH
optimization (\code{\link{maxBHHH}}) and utilities that extract
standard errors from the estimates.  Other highlights include a unified
interface for all included optimizers, tools to test user-provided analytic
derivatives, and constrained optimization.

A good starting point to learn about the usage of \pkg{maxLik} are the
included vignettes \dQuote{Introduction: what is maximum likelihood},
\dQuote{Maximum likelihood estimation with maxLik} and
\dQuote{Stochastic Gradient Ascent in maxLik}.  Another good
source is
Henningsen & Toomet (2011), an introductory paper to the package.
Use
\code{vignette(package="maxLik")} to see the available vignettes, and
\code{vignette("using-maxlik")} to read the usage vignette.

From the user's perspective, the
central function in the package is \code{\link{maxLik}}.  In its
simplest form it takes two arguments: the log-likelihood function, and
a vector of initial parameter values (see the example below).
It returns an object of class
\sQuote{maxLik} with convenient methods such as
\code{\link[=summary.maxLik]{summary}},
\code{\link[=coef.maxLik]{coef}}, and
\code{\link[=stdEr.maxLik]{stdEr}}.  It also supports a plethora
of other arguments, for instance one can supply analytic gradient and
Hessian, select the desired optimizer, and control the optimization in
different ways.

A useful utility functions in the package is
\code{\link{compareDerivatives}} that
allows one to compare the analytic and numeric derivatives for debugging
purposes.
Another useful function is \code{\link{condiNumber}} for
analyzing multicollinearity problems in the estimated models.

In the interest of providing a unified user interface, all the
optimizers are implemented as maximizers in this package.  This includes
the \code{\link{optim}}-based methods, such as \code{\link{maxBFGS}} and
\code{\link{maxSGA}}, the maximizer version of popular Stochastic
Gradient Descent.
}
\author{
Ott Toomet <otoomet@gmail.com>,
Arne Henningsen <arne.henningsen@gmail.com>,
with contributions from Spencer Graves, Yves Croissant and David
Hugh-Jones. 

Maintainer: Ott Toomet <otoomet@gmail.com>
}
\references{
  Henningsen A, Toomet O (2011). \dQuote{maxLik: A package for maximum
    likelihood estimation in R.}
  Computational Statistics, 26(3), 443-458. doi:
  \doi{10.1007/s00180-010-0217-1}.
}
\keyword{Basics|package}
\keyword{Mathematics|optimize}
\examples{
### estimate mean and variance of normal random vector

## create random numbers where mu=1, sd=2
set.seed(123)
x <- rnorm(50, 1, 2 )

## log likelihood function.
## Note: 'param' is a 2-vector c(mu, sd)
llf <- function(param) {
   mu <- param[1]
   sd <- param[2]
   llValue <- dnorm(x, mean=mu, sd=sd, log=TRUE)
   sum(llValue)
}

## Estimate it with mu=0, sd=1 as start values
ml <- maxLik(llf, start = c(mu=0, sigma=1) )
print(summary(ml))
## Estimates close to c(1,2) :-)
}
