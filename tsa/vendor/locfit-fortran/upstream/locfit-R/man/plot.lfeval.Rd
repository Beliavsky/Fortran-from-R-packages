\name{plot.lfeval}
\alias{plot.lfeval}
\title{
Plot a Locfit Evaluation Structure.
}
\usage{
\method{plot}{lfeval}(x, add=FALSE, txt=FALSE, ...)
}
\description{
  Plots the evaluation points from a \code{locfit} or \code{lfeval}
  structure, for one- or two-dimensional fits.
}
\arguments{
  \item{x}{A \code{lfeval} or \code{locfit} object}
  \item{add}{If \code{TRUE}, the points will be added to the existing plot.
    Otherwise, a new plot is created.}
  \item{txt}{If \code{TRUE}, the points are annotated with numbers in the
    order they were entered into the fit.}
  \item{...}{Additional graphical parameters.}
}

\value{
  \code{"lfeval"} object.
}

\seealso{
  \code{\link{lfeval}},
  \code{\link{locfit}},
  \code{\link{print.lfeval}}
}

\keyword{smooth}
