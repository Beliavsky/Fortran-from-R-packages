\name{lines.polynomial}
\alias{lines.polynomial}
\title{Lines Method for Polynomials}
\description{
  Add a polynomial to an existing plot usually as a line plot.
}
\usage{
\method{lines}{polynomial}(x, len = 1000, xlim = NULL, ylim = NULL, \dots)
}
\arguments{
  \item{x}{an object of class \code{"polynomial"}.}
  \item{len}{size of vector at which evaluations are to be made.}
  \item{xlim, ylim}{the range of x and y values with sensible defaults.}  
  \item{\dots}{additional arguments as for the \code{\link{lines}}
    generic.}
}
\details{
  This is a method for the generic function \code{\link{lines}}.
    
  Lines representing the given polynomial are added to an existing
  plot.  Values outside the current plot region are not shown.
}
\seealso{
  \code{\link{lines}},
  \code{\link{points}},
  \code{\link{points.polynomial}},
  \code{\link{plot}},
  \code{\link{plot.polynomial}}.
}
\examples{
plot  (poly.calc(-1:5))
lines (poly.calc( 2:4), lty = 2)
points(poly.calc(-2:6), pch = 4)
}
\keyword{symbolmath}
