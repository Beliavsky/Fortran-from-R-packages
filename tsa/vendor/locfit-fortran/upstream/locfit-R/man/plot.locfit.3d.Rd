\name{plot.locfit.3d}
\alias{plot.locfit.3d}
\title{
Plot a high-dimensional "preplot.locfit" object using trellis displays.
}
\usage{
\method{plot}{locfit.3d}(x, main="", pv, tv, type = "level", pred.lab = x$vnames,
               resp.lab=x$yname, crit = 1.96, ...)
}
\description{
This function plots cross-sections of a Locfit model (usually in three
or more dimensions) using trellis displays. It is not usually called
directly, but is invoked by \code{\link{plot.locfit}}.

The R libraries \code{lattice} and \code{grid} provide a partial
(at time of writing) implementation of trellis. Currently, this works
with one panel variable.
}
\arguments{
  \item{x}{\code{"preplot.locfit"} object.}
  \item{main}{title for the plot.}
  \item{pv}{Panel variables. These are the variables (either one or two)
    that are varied within each panel of the display.}
  \item{tv}{Trellis variables. These are varied from panel to panel
    of the display.}
  \item{type}{Type of display. When there are two panel variables,
    the choices are \code{"contour"}, \code{"level"} and
    \code{"persp"}.}
  \item{pred.lab}{label for the predictor variable.}
  \item{resp.lab}{label for the response variable.}
  \item{crit}{critical value for the confidence level.}
  \item{...}{graphical parameters passed to \code{xyplot} or
    \code{contourplot}.} 
}
\seealso{
  \code{plot.locfit},
  \code{preplot.locfit}
}

\keyword{methods}
