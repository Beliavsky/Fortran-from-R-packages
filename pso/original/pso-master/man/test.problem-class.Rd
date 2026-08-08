\name{test.problem-class}
\Rdversion{1.1}
\docType{class}
\alias{test.problem-class}

\title{Class "test.problem"}
\description{
The class contains a test problem including domain
definition and reference solution. Generally objects from the class
facilitate easy testing of PSO with various parameters.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("test.problem",
  ...)}, but the convenience constructer \code{\link{test.problem}} is
the usual approach.
}
\section{Slots}{
  \describe{
    \item{\code{name}:}{The name of the test problem. Object of class \code{"character"}.}
    \item{\code{f}:}{Function to be minimized. Object of class \code{"function"}.}
    \item{\code{grad}:}{Gradient of \code{f}. Only used with BFGS for
      the local search. Object of class \code{"function"}.}
    \item{\code{n}:}{Problem dimensionality. Object of class \code{"integer"}.}
    \item{\code{maxf}:}{Maximal number of function evaluations to use. Object of class \code{"integer"}}
    \item{\code{objective}:}{The absolute tolerance when running PSO. Object of class \code{"numeric"}.}
    \item{\code{ntest}:}{The number of tests to perform. Object of class \code{"integer"}.}
    \item{\code{lower}:}{The lower bounds. Object of class \code{"numeric"}.}
    \item{\code{upper}:}{The upper bounds. Object of class \code{"numeric"}.}
  }
}
\section{Methods}{
  \describe{
    \item{psoptim}{\code{signature(par = "test.problem", fn = "missing",
	gr = "missing",}\cr
\code{lower = "missing", upper = "missing")}: for
      running PSO on the test problem. See \code{\link{psoptim-methods}}
      for details. }
    \item{show}{\code{signature(object = "test.problem")}: descriptive
      information of the test problem. See \code{\link{show-methods}}
      for details.}
	 }
}

\seealso{
  \code{\link{test.problem}}.
}
\examples{
test.problem("rast")
test.problem("rast",10) # modified for 10 repetitions.

test.problem("para")
}
\keyword{classes}
