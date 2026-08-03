\name{summary.maxLik}
\alias{summary.maxLik}
\alias{coef.summary.maxLik}
\title{summary the Maximum-Likelihood estimation}
\description{
  Summary the Maximum-Likelihood estimation including standard errors
  and t-values.
}
\usage{
\method{summary}{maxLik}(object, eigentol=1e-12, ... )
\method{coef}{summary.maxLik}(object, \ldots)
}

\arguments{
  \item{object}{
    object of class 'maxLik', or 'summary.maxLik',
    usually a result from
    Maximum-Likelihood estimation.
  }
  \item{eigentol}{
    The standard errors are only calculated if the ratio of the smallest
    and largest eigenvalue of the Hessian matrix is less than
    \dQuote{eigentol}.  Otherwise the Hessian is treated as singular.
  }
  \item{\ldots}{currently not used.}
}
\value{
  An object of class 'summary.maxLik' with following components:
  \describe{
    \item{type}{type of maximization.}
    \item{iterations}{number of iterations.}
    \item{code}{code of success.}
    \item{message}{a short message describing the code.}
    \item{loglik}{the loglik value in the maximum.}
    \item{estimate}{numeric matrix, the first column contains the parameter
      estimates, the second the standard errors, third t-values and fourth
      corresponding probabilities.}
    \item{fixed}{logical vector, which parameters are treated as constants.}
    \item{NActivePar}{number of free parameters.}
    \item{constraints}{information about the constrained optimization.
      Passed directly further from \code{maxim}-object.  \code{NULL} if
      unconstrained maximization.
    }
  }
}
\author{Ott Toomet, Arne Henningsen}
\seealso{
  \code{\link{maxLik}} for maximum likelihood estimation,
  \code{\link{confint}} for confidence intervals, and \code{\link{tidy}}
  and \code{\link{glance}} for alternative quick summaries of the ML
  results. 
}
\examples{
## ML estimation of exponential distribution:
t <- rexp(100, 2)
loglik <- function(theta) log(theta) - theta*t
gradlik <- function(theta) 1/theta - t
hesslik <- function(theta) -100/theta^2
## Estimate with numeric gradient and hessian
a <- maxLik(loglik, start=1, control=list(printLevel=2))
summary(a)
## Estimate with analytic gradient and hessian
a <- maxLik(loglik, gradlik, hesslik, start=1, control=list(printLevel=2))
summary(a)
}
\keyword{models}
