\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,test.result,missing-method}
\title{Plot methods for test.result objects}
\description{
Graphical methods for plotting test results.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "test.result", y = "missing")}}{
Produces a plot of the success rate versus the number of function
evaluations for the test result provided as \code{x}. Any additional
arguments to the method will be passed on to \code{\link{plot}}.
}
}}
\seealso{
\code{\link{plot}}, \code{\linkS4class{test.result}}.
}
\keyword{methods}
