\name{plot.gpd}
\alias{plot.gpd}
\title{Plot Fitted GPD Model}
\description{
  The plot method plot.gpd provides four different plots 
  for assessing fitted GPD model. The user selects the plot
  type from a menu. See the examples below. 
}
\usage{
\method{plot}{gpd}(x, optlog = NA, extend = 1.5, labels = TRUE, \dots)
}
\arguments{
  \item{x}{a \code{gpd} object}
  \item{optlog}{optional argument for plots 1 and 2 giving a
    particular choice of logarithmic axes: \code{"x"} x-axis
    only; \code{"y"} y-axis only; \code{"xy"} both axes;
    \code{""} neither axis.}
  \item{extend}{optional argument for plots 1 and 2 expressing how
    far x-axis should extend as a multiple of the largest data value.
    This argument must take values greater than 1 and is useful for
    showing estimated quantiles beyond data.}
  \item{labels}{optional argument for plots 1 and 2  specifying whether
    or not axes should be labelled}
  \item{\dots}{other graphics parameters}
}
\value{
  If plot 1 or 2 is selected as the final plot, a list object containing
  details of the plot is returned invisibly.
  This object should be used as the first argument of \code{gpd.q}
  or \code{gpd.sfall} to add quantile estimates or expected shortfall
  estimates to the plot.
}
\seealso{\code{\link{gpd}}, \code{\link{quant}}, \code{\link{shape}}}

\examples{
data(danish)
out <- gpd(danish, 10) 
\dontrun{plot(out)} 

\dontrun{Make a plot selection (or 0 to exit):}
\dontrun{1: plot: Excess Distribution} 
\dontrun{2: plot: Tail of Underlying Distribution} 
\dontrun{3: plot: Scatterplot of Residuals} 
\dontrun{4: plot: QQplot of Residuals} 
}
\keyword{hplot}
