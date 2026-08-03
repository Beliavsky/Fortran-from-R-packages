\name{storedValues}
\alias{storedValues}
\alias{storedValues.maxim}
\alias{storedParameters}
\alias{storedParameters.maxim}
\title{Return the stored values of optimization} 
\description{
  Retrieve the objective function value for each iteration if stored
  during the optimization.
}
\usage{
storedValues(x, \dots)
\method{storedValues}{maxim}(x, \dots)
storedParameters(x, \dots)
\method{storedParameters}{maxim}(x, \dots)
}
\arguments{
  \item{x}{a result of maximization, created
    by \code{\link{maxLik}}, \code{\link{maxSGA}} or another optimizer.}
  \item{\dots}{further arguments for other methods}
}
\details{
  These is a generic method.
  If asked by control parameter \code{storeValues=TRUE} or \code{storeParameters=TRUE}, certain
  optimization methods store the objective function value and the
  parameter value at each epoch.  These methods
  retrieves the stored values.
}
\value{
  \itemize{
    \item \code{storedValues}: a numeric vector, one value for each
    iteration
    \item \code{storedParameters}: a numeric matrix with rows
    corresponding to the iterations and columns to the parameter
    components.
  }
  In both cases, the first value stored corresponds to the initial
  parameter. 
}
\author{Ott Toomet}
\seealso{\code{\link{maxSGA}},
  \code{\link{maxControl}}
  }
\examples{
## Estimate the exponential distribution parameter
t <- rexp(100, 2)
loglik <- function(theta, index) sum(log(theta) - theta*t[index])
## Estimate with numeric gradient and numeric Hessian
a <- maxSGA(loglik, start=1,
            control=list(storeValues=TRUE, storeParameters=TRUE, iterlim=10),
            nObs=100)
storedValues(a)
storedParameters(a)
}
\keyword{methods}
