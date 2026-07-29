\name{plot.gev}
\alias{plot.gev}
\title{Plot Fitted GEV Model}
\description{
  The plot method plot.gev provides two different residual plots 
  for assessing fitted GEV model. The user selects the plot type 
  from a menu. See the examples below. 
}
\usage{
\method{plot}{gev}(x, \dots)
}
\arguments{
  \item{x}{a \code{gev} object}
  \item{\dots}{other graphics parameters}
}
\details{
  Data are converted to unit exponentially distributed residuals
  under null hypothesis that GEV fits. Two diagnostics for iid
  exponential data are offered. 
}
\seealso{\code{\link{gev}}, \code{\link{qplot}}}

\examples{
data(bmw)
out <- gev(bmw, 100) 
\dontrun{plot(out)}
 
\dontrun{Make a plot selection (or 0 to exit):}
\dontrun{1: plot: Scatterplot of Residuals}
\dontrun{2: plot: QQplot of Residuals}
}
\keyword{hplot}
