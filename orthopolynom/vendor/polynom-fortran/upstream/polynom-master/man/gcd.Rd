\name{gcd}
\alias{GCD}
\alias{GCD.polylist}
\alias{GCD.polynomial}
\alias{LCM}
\alias{LCM.polylist}
\alias{LCM.polynomial}
\title{GCD and LCM for Polynomials}
\description{
  Compute the greatest common divisor (GCD) and least common multiple
  (LCM) of a collection of polynomials and polylists.
}
\usage{
\method{GCD}{polylist}(\dots)
\method{GCD}{polynomial}(\dots)
\method{LCM}{polylist}(\dots)
\method{LCM}{polynomial}(\dots)
}
\arguments{
  \item{\dots}{a list of objects of class \code{\link{polynomial}} or
    \code{\link{polylist}}.}
}
\examples{
pl <- polylist(poly.from.roots(-1),
               poly.from.roots(c(-1, -1)),
               poly.from.roots(1))
GCD(pl)
GCD(pl[-3])
LCM(pl)
LCM(pl, pl, pl[[2]])
}
\keyword{symbolmath}
