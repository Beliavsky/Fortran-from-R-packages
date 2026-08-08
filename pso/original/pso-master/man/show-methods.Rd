\name{show-methods}
\docType{methods}
\alias{show-methods}
\alias{show,test.problem-method}
\alias{show,test.result-method}
\title{ Methods for Function show }
\description{
Displays descriptive information of the object provided as argument.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "test.problem")}}{
Provide information on test problem. This includes: problem name,
dimension, objective value, maximal number of function evaluations, and
the number of test repetitions to perform.
}

\item{\code{signature(object = "test.result")}}{
Provide summary statistics for the test. This includes information on
the mean, s.d., min and max obtained for the value over all
conducted repetitions
as well as the overall success rate (percentage of test runs for
which the target objective was reached) and a measure of efficiency (the
area under the success-rate curve normalized to the maximal area
possible). Additionally displays the timing information for the test conducted.
}
}}
\keyword{methods}
