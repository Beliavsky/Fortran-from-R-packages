%% $Id: splineKnots.Rd,v 1.2 1999/10/29 22:48:32 bates Exp $
\name{splineKnots}
\title{Knot Vector from a Spline}
\usage{
splineKnots(object)
}
\alias{splineKnots}
\alias{splineKnots.spline}
\arguments{
 \item{object}{an object that inherits from class \code{spline}.}
}
\description{
  Return the knot vector corresponding to a spline object.
}
\value{
  A non-decreasing numeric vector of knot positions.
}
\author{Douglas Bates and Bill Venables}
\examples{
library( splines )
data( women )
ispl <- interpSpline( weight ~ height, women )
splineKnots( ispl )
}
\keyword{ models }
