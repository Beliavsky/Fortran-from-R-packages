\name{plot.gpdbiv}
\alias{plot.gpdbiv}
\title{Plot Fitted Bivariate GPD Model}
\description{
  Provides a number of plots summarising a bivariate GPD model
  fitted using the bivariate POT method. See the examples below. 
}
\usage{
\method{plot}{gpdbiv}(x, extend = 1.1, n.contours = 15, \dots)
}
\arguments{
  \item{x}{a \code{gpdbiv} object}
  \item{extend}{optional argument expressing how far x-axis should
    extend as a multiple of the largest data value.}
  \item{n.contours}{number of contours in bivariate contour plots}
  \item{\dots}{other graphics parameters}
}
\details{
  Option 1 plots the threshold exceedance data; option 2 plots contours
  of the fitted bivariate distribution function in the joint upper tail
  (above both thresholds); option 3 plots corresponding contours of the
  fitted joint survival function; plots 4 and 5 show the fitted tails
  of the marginal distributions. 
}
\seealso{\code{\link{gpd}}, \code{\link{gpdbiv}},
  \code{\link{tailplot}}, \code{\link{interpret.gpdbiv}},
  \code{\link{plot.gpd}}}

\examples{
data(bmw) ; data(siemens)
out <- gpdbiv(-bmw, -siemens, ne1 = 100, ne2 = 100) 
\dontrun{plot(out)}
 
\dontrun{Make a plot selection (or 0 to exit): }
\dontrun{1: plot: Exceedance data}
\dontrun{2: plot: Contours of Bivariate Distribution Function}
\dontrun{3: plot: Contours of Bivariate Survival Function} 
\dontrun{4: plot: Tail of Marginal 1} 
\dontrun{5: plot: Tail of Marginal 2} 
}
\keyword{hplot}
