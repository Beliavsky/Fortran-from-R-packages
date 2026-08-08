\name{test.result-class}
\Rdversion{1.1}
\docType{class}
\alias{test.result-class}

\title{Class "test.result"}
\description{
  A container class with results from executing a (repetition of) test
  problem(s).
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("test.result",
  ...)}, but the object is normally provided as the result of executing
\code{psoptim} on an object of class \code{"test.problem"}.
}
\section{Slots}{
  \describe{
    \item{\code{problem}:}{Object of class \code{"test.problem"}.}
    \item{\code{result}:}{A list with each of the results from
      repetitive invocation of \code{\link{psoptim}} on \code{problem}. }
    \item{\code{time}:}{The overall time taken for executing the test.
      Object of class \code{"numeric"}. }
  }
}
\section{Methods}{
  \describe{
    \item{getSuccessRate}{\code{signature(object = "test.result")}:
      internal method used to calculate the success rate for a series of
      test results. See \code{\link{getSuccessRate-methods}}
      for details.}
    \item{lines}{\code{signature(x = "test.result")}: add lines with the
    test result to an existing plot. See \code{\link{lines-methods}}
      for details.}
    \item{plot}{\code{signature(x = "test.result", y = "missing")}: plot
    the test result. See \code{\link{plot-methods}}
      for details.}
    \item{points}{\code{signature(x = "test.result")}: add points with
      the test result to an existing plot. See \code{\link{points-methods}}
      for details.}
    \item{show}{\code{signature(object = "test.result")}: summary
      statistics of the test. See \code{\link{show-methods}}
      for details.}
	 }
}
\examples{
showClass("test.result")

set.seed(1)
t <- test.problem("rastrigin",10)
o <- psoptim(t)
show(o)

\dontrun{
plot(o)
}
}
\keyword{classes}
