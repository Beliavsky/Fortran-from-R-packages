\name{preplot.locfit.raw}
\alias{preplot.locfit.raw}
\title{
Prediction from a Locfit object.
}
\usage{
\method{preplot}{locfit.raw}(object, newdata, where, what, band, ...)
}
\arguments{
\item{object}{
Fitted object from \code{\link{locfit}()}.
}
\item{newdata}{
New data points.
}
\item{where}{
Type of data provided in \code{newdata}.
}
\item{what}{
What to compute predicted values of.
}
\item{band}{
  Compute standard errors for the fit and include confidence bands on
  the returned object.}
\item{...}{Arguments passed to and from other methods.}
}
\value{
A list containing raw output from the internal prediction routines.
}
\description{
  \code{preplot.locfit.raw} is an internal function used by
  \code{\link{predict.locfit}} and \code{\link{preplot.locfit}}.
  It should not normally be called directly.
}
\seealso{
\code{\link{locfit}}, \code{\link{predict.locfit}}, \code{\link{preplot.locfit}}.
}
\keyword{methods}
% Converted by Sd2Rd version 0.2-a5.
