\name{AIC.maxLik}
\alias{AIC.maxLik}
\alias{coef.maxim}
\alias{coef.maxLik}
\alias{stdEr.maxLik}

\title{Methods for the various standard functions}
\description{
  These are methods for the maxLik related objects.  See also the
  documentation for the
  corresponding generic functions
}
\usage{
\method{AIC}{maxLik}(object, \dots, k=2)
\method{coef}{maxim}(object, \dots)
\method{coef}{maxLik}(object, \dots)
\method{stdEr}{maxLik}(x, eigentol=1e-12, \dots)
}
\arguments{
  \item{object}{a \sQuote{maxLik} object (\code{coef} can also handle
    \sQuote{maxim} objects)}
  \item{k}{numeric, the penalty per parameter to be used; the default
    \sQuote{k = 2} is the classical AIC.}
  \item{x}{a \sQuote{maxLik} object}
  \item{eigentol}{
    The standard errors are only calculated if the ratio of the smallest
    and largest eigenvalue of the Hessian matrix is less than
    \dQuote{eigentol}.  Otherwise the Hessian is treated as singular.
  }
  \item{\dots}{other arguments for methods}
}
\details{
  \describe{
    \item{AIC}{calculates Akaike's Information Criterion (and other
      information criteria).}
    \item{coef}{extracts the estimated parameters (model's
      coefficients).}
    \item{stdEr}{extracts standard errors (using the Hessian matrix).
    }
  }
}
\examples{
## estimate mean and variance of normal random vector
set.seed(123)
x <- rnorm(50, 1, 2)

## log likelihood function.
## Note: 'param' is a vector
llf <- function( param ) {
   mu <- param[ 1 ]
   sigma <- param[ 2 ]
   return(sum(dnorm(x, mean=mu, sd=sigma, log=TRUE)))
}

## Estimate it.  Take standard normal as start values
ml <- maxLik(llf, start = c(mu=0, sigma=1) )

coef(ml)
stdEr(ml)
AIC(ml)
}
\keyword{methods}
