\name{meplot}
\alias{meplot}
\title{Sample Mean Excess Plot}
\description{
  Plots sample mean excesses over increasing thresholds.  
}
\usage{
meplot(data, omit = 3, labels = TRUE, \dots)
}
\arguments{
  \item{data}{data vector}
  \item{omit}{number of upper plotting points to be omitted}
  \item{labels}{whether or not axes are to be labelled}
  \item{\dots}{other graphics parameters}
}
\details{
  An upward trend in plot shows heavy-tailed behaviour. In particular,
  a straight line with positive gradient above some threshold is a sign 
  of Pareto behaviour in tail. A downward trend shows thin-tailed 
  behaviour whereas a line with zero gradient shows an exponential 
  tail. Because upper plotting points are the average of a handful of
  extreme excesses, these may be omitted for a prettier plot. 
}
\seealso{\code{\link{gpd}}, \code{\link{qplot}}}

\examples{
\dontrun{data(danish)}
\dontrun{meplot(danish)} 
# Sample mean excess plot of heavy-tailed Danish fire insurance data 
}
\keyword{hplot}
