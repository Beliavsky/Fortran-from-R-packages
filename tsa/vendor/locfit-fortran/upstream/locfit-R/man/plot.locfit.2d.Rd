\name{plot.locfit.2d}
\alias{plot.locfit.2d}
\title{
Plot a two-dimensional "preplot.locfit" object.
}
\usage{
\method{plot}{locfit.2d}(x, type="contour", main, xlab, ylab, zlab=x$yname, ...)
}
\arguments{
  \item{x}{Two dimensional \code{preplot.locfit} object.}
  \item{type}{one of \code{"contour"}, \code{"persp"}, or \code{"image"}.}
  \item{main}{title for the plot.}
  \item{xlab, ylab}{text labels for the x- and y-axes.}
  \item{zlab}{if \code{type="persp"}, the label for the z-axis.}
  \item{...}{Additional arguments to the \code{contour}, \code{persp} or
    \code{image} functions.}
}

\description{
This function is not usually called directly. It will be called automatically
when plotting one-dimensional \code{locfit} or \code{preplot.locfit}
objects.
}

\seealso{
\code{\link{locfit}}, \code{\link{plot.locfit}}, \code{\link{preplot.locfit}}
}
\keyword{methods}
