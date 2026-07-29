\name{qplot}
\alias{qplot}
\title{Exploratory QQplot for Extreme Value Analysis}
\description{
  Creates a QQplot for threshold data against the exponential 
  distribution or the generalized Pareto distribution. 
}
\usage{
qplot(data, xi = 0, trim = NA, threshold = NA, line = TRUE,
    labels = TRUE, \dots)
}
\arguments{
  \item{data}{data vector}
  \item{xi}{the xi value of a generalized Pareto distribution}
  \item{trim}{value at which data are to be right-truncated}
  \item{threshold}{value at which data are to be left-truncated}
  \item{line}{whether or not a straight line is to be added}
  \item{labels}{whether or not the axes are to be labelled}
  \item{\dots}{other graphics parameters}
}
\details{
  If xi is zero the reference distribution is the exponential; if 
  xi is non-zero the reference distribution is the generalized 
  Pareto with that value of xi. 
  In the case of the exponential, the plot is interpreted as follows. 
  Concave departures from a straight line are a sign of heavy-tailed  
  behaviour. Convex departures show thin-tailed behaviour. 
}
\seealso{\code{\link{gpd}}, \code{\link{meplot}}}

\examples{
\dontrun{data(danish)}
\dontrun{qplot(danish)} 
# QQplot of heavy-tailed Danish fire insurance data 
}
\keyword{hplot}
