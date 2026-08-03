\name{activePar}
\alias{activePar}
\alias{activePar.default}
\title{free parameters under maximization}
\description{
  Return a logical vector, indicating which parameters were free under
  maximization, as opposed to the fixed parameters that are treated as
  constants.  See argument \dQuote{fixed} for \code{\link{maxNR}}.
}
\usage{
activePar(x, \dots)
\method{activePar}{default}(x, \dots)
}
\arguments{
  \item{x}{object, created by a maximization routine, such as
    \code{\link{maxNR}} or \code{\link{maxLik}}, or derived from a
    maximization object.
  }
  \item{\dots}{further arguments for methods}
}
\details{
  Several optimization routines allow the user to fix some parameter
  values (or do it automatically in some cases).  For gradient or
  Hessian based inference one has to know which parameters carry
  optimization-related information.
}
\value{
  A logical vector, indicating whether the parameters were free to
  change during optimization algorithm.
}
\author{Ott Toomet}
\seealso{\code{\link{maxNR}}, \code{\link[miscTools]{nObs}}}
\examples{
## a two-dimensional exponential hat
f <- function(a) exp(-a[1]^2 - a[2]^2)

## maximize wrt. both parameters 
free <- maxNR(f, start=1:2) 
summary(free)  # results should be close to (0,0)
activePar(free)

## keep the first parameter constant
cons <- maxNR(f, start=1:2, fixed=c(TRUE,FALSE))
summary(cons) # result should be around (1,0)
activePar(cons)
}
\keyword{methods}
\keyword{optimize}
