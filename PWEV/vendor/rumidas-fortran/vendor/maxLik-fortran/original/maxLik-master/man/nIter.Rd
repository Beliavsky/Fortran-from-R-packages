\name{nIter}
\alias{nIter}
\alias{nIter.default}
\title{Return number of iterations for iterative models} 
\description{
  Returns the number of iterations for iterative models.  The default
  method assumes presence of a component \code{iterations} in \code{x}.
}
\usage{
nIter(x, \dots)
\method{nIter}{default}(x, \dots)
}
\arguments{
  \item{x}{a statistical model, or a result of maximisation, created
    by \code{\link{maxLik}}, \code{\link{maxNR}} or another optimizer.}
  \item{\dots}{further arguments for methods}
}
\details{
  This is a generic function.  The default method returns the component
  \code{x$iterations}.
}
\value{
  numeric, number of iterations.  Note that \sQuote{iteration} may mean
  different things for different optimizers.
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
nIter(a)
}
\keyword{methods}
