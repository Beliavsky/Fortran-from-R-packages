\name{points-methods}
\docType{methods}
\alias{points-methods}
\alias{points,test.result-method}
\title{Methods for Function points}
\description{
Graphical methods for adding points to existing plots.
}
\section{Methods}{
\describe{

\item{\code{signature(x = "test.result")}}{
Add points with the success rate versus the number of function
evaluations for the test resulted provided as \code{x} to the current
plot. Any additional 
arguments to the method will be passed on to
\code{\link{points}}. Typically this method is used to add new test
results to an existing plot.
}
}}
\seealso{
\code{\link{points}}, \code{\linkS4class{test.result}}.
}
\keyword{methods}

