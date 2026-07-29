\name{emplot}
\alias{emplot}
\title{Plot of Empirical Distribution Function}
\description{
  Plots empirical distribution function of a sample. 
}
\usage{
emplot(data, alog = "x", labels = TRUE, \dots)
}
\arguments{
  \item{data}{data vector}
  \item{alog}{whether axes are to be logged: \code{"x"} x-axis
    only; \code{"y"} y-axis only; \code{"xy"} both axes;
    \code{""} neither axis.}
  \item{labels}{whether or not axes should be labelled}
  \item{\dots}{other graphics parameters}
}
\details{
  This is a simple explanatory function. A straight line on the  
  double log scale indicates Pareto tail behaviour. 
}
\seealso{\code{\link{qplot}}, \code{\link{meplot}}}

\examples{
\dontrun{data(danish)}
\dontrun{emplot(danish)} 
# Danish fire insurance data show Pareto tail behaviour 
}
\keyword{hplot}
