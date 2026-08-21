\name{points.locfit}
\alias{points.locfit}
\title{
Add `locfit' points to existing plot
}
\usage{
\method{points}{locfit}(x, tr, ...)
}
\arguments{
\item{x}{
\code{"locfit"} object. Should be a model with one predictor.
}
\item{tr}{
Back transformation.
}
\item{...}{
Other arguments to the default \code{\link{points}} function.
}}
\description{This function shows the points at which the local fit
  was computed directly, rather than being interpolated. This can
  be useful if one is unsure of the validity of interpolation.}
\seealso{
\code{\link{locfit}}, \code{\link{plot.locfit}}, \code{\link{points}}
}
\keyword{smooth}
% Converted by Sd2Rd version 0.2-a5.
