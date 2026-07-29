\name{step_plot}
\alias{step_plot}
\alias{edf_plot}
\alias{eqf_plot}
\title{Plot of Step Functions, Empirical Distribution and Quantile Functions}
\description{
  Plotting step functions, empirical distribution functions and
  empirical quantile functions.
}
\usage{
step_plot(x, y, y0 = NA, x0 = NA, x1 = NA, method = c("edf", "eqf"), log = "",
          verticals = NA, do.points = NA, add = FALSE,
          col = par("col"), main = "", xlab = "x", ylab = "Function value at x",
          plot.args = NULL, segments.args = NULL, points.args = NULL)
edf_plot(x, y0 = 0, x0 = NA, x1 = NA, log = "",
         verticals = NA, do.points = NA, col = par("col"),
         main = "", xlab = "x", ylab = "Distribution function at x", ...)
eqf_plot(x, y0 = NA, x0 = 0, x1 = 1, log = "",
         verticals = NA, do.points = NA, col = par("col"),
         main = "", xlab = "x", ylab = "Quantile function at x", ...)
}
\arguments{
  \item{x}{
    \describe{
      \item{\code{step_plot()}:}{\code{\link{numeric}} vector of
	x-values.}
      \item{\code{edf_plot()}:}{\code{\link{numeric}} vector or a
	\code{\link{list}} of \code{\link{numeric}} vectors; if a list, each
	element corresponds to the x-values of an empirical distribution function.}
      \item{\code{eqf_plot()}:}{similar to \code{edf_plot()}.}
    }
  }
  \item{y}{y-values corresponding to \code{x}.}
  \item{y0}{y-value of the graph extending to the left of the
    first x-value.}
  \item{x0}{smallest x-value.}
  \item{x1}{largest x-value.}
  \item{method}{\code{\link{character}} string indicating the type of
    method to be used (\code{"edf"} for empricial distribution function
    types of plots and \code{"eqf"} for empirical quantile function types).}
  \item{log}{\code{\link{character}} indicating whether
     a logarithmic x-axis is used.}
  \item{verticals}{\code{\link{logical}} indicating whether to plot
    vertical lines (defaults to \code{\link{TRUE}} if and only if
    there are 100 or more data points).}
  \item{do.points}{\code{\link{logical}} (vector) indicating whether points are
    to be plotted (defaults to \code{\link{TRUE}} if and only if there
    are less than 100 data points).}
  \item{add}{\code{\link{logical}} indicating whether the current plot
    is added to the last one.}
  \item{col}{color (for \code{edf_plot()} this can be a
    \code{\link{vector}}).}
  \item{main}{title.}
  \item{xlab}{x-axis label.}
  \item{ylab}{y-axis label.}
  \item{plot.args}{\code{\link{list}} of additional arguments passed to
    the underlying \code{\link{plot}()}.}
  \item{segments.args}{\code{\link{list}} of additional arguments passed to
    the underlying \code{\link{segments}()}.}
  \item{points.args}{\code{\link{list}} of additional arguments passed to
    the underlying \code{\link{points}()}.}
  \item{\dots}{additional arguments passed to the underlying
    \code{step_plot()}.}
}
\value{
  Nothing (plot by side-effect).
}
\author{Marius Hofert}
\examples{
x <- c(5, 2, 4, 2, 3, 2, 2, 2, 1, 2) # example data
edf_plot(x) # empirical distribution function (edf)
edf_plot(x, log = "x")
edf_plot(x, verticals = TRUE)
edf_plot(x, do.points = FALSE)
cols <- c("black", "royalblue3")
edf_plot(list(x, x+2), col = cols) # edf with shifted edf
edf_plot(list(x, x+2), col = cols, x0 = 0.5, x1 = 7.5)
edf_plot(list(x, x+2), col = cols, x0 = 0.5, x1 = 7.5, verticals = TRUE)
eqf_plot(x) # empirical quantile function
eqf_plot(x, verticals = TRUE)
}
\keyword{hplot}