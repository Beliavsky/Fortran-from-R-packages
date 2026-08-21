\name{plot.locfit}
\alias{plot.locfit}
\title{
Plot an object of class locfit.
}
\usage{
\method{plot}{locfit}(x, xlim, pv, tv, m, mtv=6, band="none", tr=NULL,
  what = "coef", get.data=FALSE, f3d=(d == 2) && (length(tv) > 0), ...)
}
\arguments{
\item{x}{
locfit object.
}
\item{xlim}{
Plotting limits. Eg. \code{xlim=c(0,0,1,1)} plots over the unit square in
two dimensions. Default is bounding box of the data.
}
\item{pv}{
Panel variables, to be varied within each panel of a plot.  May be
        specified as a character vector, or variable numbers.  There must be
        one or two panel variables; default is all variables in one or two
        dimensions; Variable 1 in three or more dimensions.
May by specified using either variable numbers or names.
}
\item{tv}{
Trellis variables, to be varied from panel to panel of the plot.
}
\item{m}{
Controls the plot resolution (within panels, for trellis displays).
Default is 100 points in one dimension; 40 points (per dimension)
in two or more dimensions.
}
\item{mtv}{
Number of points for trellis variables; default 6.
}
\item{band}{
  Type of confidence bands to add to the plot. Default is \code{"none"}. Other
  choices include \code{"global"} for bands using a global variance estimate;
  \code{"local"} for bands using a local variance estimate and \code{"pred"}
  for prediction bands (at present, using a global variance estimate).
  To obtain the global variance estimate for a fit, use \code{\link{rv}}.
  This can be changed with \code{\link{rv<-}}. Confidence bands, by default,
  are 95\%, based on normal approximations and neglecting bias.
  To change the critical value or confidence level, or to obtain
  simultaneous instead of pointwise confidence,
  the critical value stored on the fit must be changed. See the
  \code{\link{kappa0}} and \code{\link{crit}} functions.
}
\item{tr}{
Transformation function to use for plotting. Default is the inverse
link function, or the identity function if derivatives are requested.
}
\item{what}{
What to plot. See \code{\link{predict.locfit}}.
}
\item{get.data}{
If \code{TRUE}, original data is added to the plot. Default: \code{FALSE}.
}
\item{f3d}{
Force the \code{locfit.3d} class on the prediction object, thereby generating
a trellis style plot. Default: \code{FALSE}, unless a \code{tv} argument is'
provided. Not available in R.
}
\item{...}{
Other arguments to \code{plot.locfit.1d}, \code{plot.locfit.2d} or
\code{plot.locfit.3d} as appropriate.
}}
\description{
The \code{plot.locfit} function generates grids of ploting points, followed
by a call to \code{\link{preplot.locfit}}. The returned object is then
passed to \code{\link{plot.locfit.1d}}, \code{\link{plot.locfit.2d}} or
\code{\link{plot.locfit.3d}} as appropriate.
}
\examples{
x <- rnorm(100)
y <- dnorm(x) + rnorm(100) / 5
plot(locfit(y~x), band="global")
x <- cbind(rnorm(100), rnorm(100))
plot(locfit(~x), type="persp")
}
\seealso{
\code{\link{locfit}}, \code{\link{plot.locfit.1d}},
\code{\link{plot.locfit.2d}}, \code{\link{plot.locfit.3d}},
\code{\link{lines.locfit}}, \code{\link{predict.locfit}},
\code{\link{preplot.locfit}}
}
\keyword{methods}
