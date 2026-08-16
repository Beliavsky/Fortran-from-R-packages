\name{points.polynomial}
\alias{points.polynomial}
\title{Points Method for Polynomials}
\description{
  Add a polynomial to an existing plot usually as a point plot.
}
\usage{
\method{points}{polynomial}(x, length = 100, \dots)
}
\arguments{
  \item{x}{an object of class \code{"polynomial"}.}
  \item{length}{size of x vector at which evaluations are to be made.}
  \item{\dots}{additional arguments as for the points generic.}
}
\details{
  This is a method for the generic function \code{\link{points}}.

  Points representing the given polynomial are added to an existing
  plot.  Values outside the current plot region are not shown.
}
\seealso{
  \code{\link{plot}},
  \code{\link{lines}},
  \code{\link{points}},
  \code{\link{plot.polynomial}},
  \code{\link{lines.polynomial}}.
}
\examples{
plot(poly.calc(-1:5))
lines(poly.calc(2:4), lty=2)
points(poly.calc(-2:6), pch=4)
}
\keyword{symbolmath}
