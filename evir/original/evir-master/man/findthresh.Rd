\name{findthresh}
\alias{findthresh}
\title{Find Threshold}
\description{
  Finds a threshold so that a given number of
  extremes lie above.  
}
\usage{
findthresh(data, ne)
}
\arguments{
  \item{data}{data vector}
  \item{ne}{vector giving number of extremes above the threshold}
}
\details{
  When the data are tied a threshold is found so that at least 
  the specified number of extremes lie above.
}
\value{
  Vector of suitable thresholds. 
}
\seealso{\code{\link{hill}}, \code{\link{gpd}},
  \code{\link{pot}}}

\examples{
# Find threshold giving (at least) fifty exceedances for Danish data
data(danish)
findthresh(danish, 50) 
}
\keyword{manip}
