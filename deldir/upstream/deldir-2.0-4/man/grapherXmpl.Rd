\name{grapherXmpl}
\alias{grapherXmpl}
\docType{data}
\title{
    grapherXmpl
}
\description{
    A data set taken from an example in the grapherator package.
    This data set demonstrates handling a data set with duplicated
    points.
}
\usage{grapherXmpl}
\format{
  A data frame with 250 observations on the following 2 variables.
  \describe{
    \item{\code{x}}{a numeric vector}
    \item{\code{y}}{a numeric vector}
  }
}
\details{
    There are 25 duplicated points, so the net number of
    observations is 225.  These data constitute a structure (named
    \code{coordinates}) generated internally in the function
    \code{addEdgesDelaunay}.  The call is to be found in the
    examples in the help file for the \code{plot.grapherator()}
    in the \code{grapherator} package.  The relevant example
    initially threw an error, revealing a bug in \code{deldir()}
    that was triggered when there were duplicated points in the data.
}
\source{
    The \code{grapherator} package,
    \url{https://CRAN.R-project.org/package=grapherator}
}
\examples{
dgX <- deldir(grapherXmpl) # Now works!!!`
}
\keyword{datasets}
