\name{logLik.maxLik}
\alias{logLik.maxLik}
\alias{logLik.summary.maxLik}
\title{Return the log likelihood value}
\description{
  Return the log likelihood value of objects of class \code{maxLik}
      and \code{summary.maxLik}.
}
\usage{
\method{logLik}{maxLik}( object, \dots )
\method{logLik}{summary.maxLik}( object, \dots )
}

\arguments{
  \item{object}{object of class \code{maxLik} or \code{summary.maxLik},
      usually a model estimated with Maximum Likelihood}
  \item{...}{additional arguments to methods}
}
\value{
  A scalar numeric, log likelihood of the estimated model.  It has
  attribute \dQuote{df}, number of free parameters.
}
\author{
   Arne Henningsen,
   Ott Toomet
}
\seealso{\code{\link{maxLik}}}
\examples{
## ML estimation of exponential duration model:
t <- rexp(100, 2)
loglik <- function(theta) log(theta) - theta*t
gradlik <- function(theta) 1/theta - t
hesslik <- function(theta) -100/theta^2
## Estimate with analytic gradient and hessian
a <- maxLik(loglik, gradlik, hesslik, start=1)
## print log likelihood value
logLik( a )
## print log likelihood value of summary object
b <- summary( a )
logLik( b )
}
\keyword{methods}
