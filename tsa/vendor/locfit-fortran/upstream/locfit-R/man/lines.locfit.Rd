\name{lines.locfit}
\alias{lines.locfit}
\alias{llines.locfit}
\title{
Add locfit line to existing plot
}
\usage{
\method{lines}{locfit}(x, m=100, tr=x$trans, \dots)
\method{llines}{locfit}(x, m=100, tr=x$trans, \dots)
}
\description{
  Adds a Locfit line to an existing plot.  \code{llines} is for use
  within a panel function for Lattice.
}
\arguments{
\item{x}{\code{locfit} object. Should be a model with one predictor.}
\item{m}{Number of points to evaluate the line at.}
\item{tr}{Transformation function to use for plotting.  Default is the
  inverse link function, or the identity function if derivatives are
  required.} 
\item{...}{Other arguments to the default \code{\link{lines}} function.}
}
\seealso{
\code{\link{locfit}}, \code{\link{plot.locfit}}, \code{\link{lines}}
}
\keyword{smooth}
