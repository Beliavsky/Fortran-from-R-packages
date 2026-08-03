\name{maxValue}
\alias{maxValue}
\alias{maxValue.maxim}
\title{Function value at maximum} 
\description{
  Returns the function value at (estimated) maximum.
}
\usage{
maxValue(x, ...)
\method{maxValue}{maxim}(x, \dots)
}
\arguments{
  \item{x}{a statistical model, or a result of maximisation, created
    by \code{\link{maxLik}}, \code{\link{maxNR}} or another optimizer.}
  \item{\dots}{further arguments for other methods}
}
\value{
  numeric, the value of the objective function at maximum.  In general,
  it is the last calculated value in case the process did not converge.
}
\author{Ott Toomet}
\seealso{\code{\link{maxLik}}, \code{\link{maxNR}}
  }
\examples{
## Estimate the exponential distribution parameter:
t <- rexp(100, 2)
loglik <- function(theta) sum(log(theta) - theta*t)
## Estimate with numeric gradient and numeric Hessian
a <- maxNR(loglik, start=1)
maxValue(a)
}
\keyword{methods}
