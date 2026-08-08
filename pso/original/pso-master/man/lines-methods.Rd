\name{lines-methods}
\docType{methods}
\alias{lines-methods}
\alias{lines,test.result-method}
\title{Methods for Function lines}
\description{
Graphical methods for adding line segments to existing plots.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "test.result")}}{
Add lines of the success rate versus the number of function
evaluations for the test resulted provided as \code{x} to the current
plot. Any additional 
arguments to the method will be passed on to
\code{\link{lines}}. Typically this method is used to add new test
results to an existing plot.
}
}}
\seealso{
\code{\link{lines}}, \code{\linkS4class{test.result}}.
}
\keyword{methods}
