\name{getSuccessRate-methods}
\docType{methods}
\alias{getSuccessRate}
\alias{getSuccessRate-methods}
\alias{getSuccessRate,test.result-method}
\title{ Methods for Function getSuccessRate }
\description{
  Provides the success rate as the result of conducting a test. Only
  implemented method is for objects of class \code{"test.result"}
}
\section{Methods}{
Calculates the success rate from the number of successful tests
conducted as a function of the number of function evaluations used.
\describe{
\item{\code{signature(object = "test.result")}}{
  This method is used internally by the graphical functions.
  Returns a list with components:

  \describe{
    \item{feval:}{The number of function evaluations.}
    \item{rate:}{The corresponding success rate (between 0 and 1).}}
}
}}
\seealso{
  \code{\linkS4class{test.result}}.
}
\keyword{methods}
\keyword{manip}
