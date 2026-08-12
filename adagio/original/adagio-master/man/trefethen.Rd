\name{SIAM test functions}
\alias{fnTrefethen}
\alias{fnWagon}
\title{
  Trefethen and Wagon Test Functions
}
\description{
  Test functions for global optimization posed for the SIAM 100-digit
  challenge in 2002 by Nick Trefethen, Oxford University, UK.
}
\usage{
fnTrefethen(p2)
fnWagon(p3)
}
\arguments{
  \item{p2}{Numerical vector of length 2.}
  \item{p3}{Numerical vector of length 3.}
}
\details{
  These are highly nonlinear and oscillating functions in two and three
  dimensions with thousands of local mimima inside the unit square resp.
  cube (i.e., [-1, 1] x [-1, 1] or [-1, 1] x [-1, 1] x [-1, 1]).
}
\value{
  Function value is a single real number.
}
\references{
  F. Bornemann, D. Laurie, S. Wagon, and J. Waldvogel (2004). The
  SIAM 100-Digit Challenge: A Study in High-Accuracy Numerical Computing.
  Society for Industrial and Applied Mathematics.
}
\author{
  HwB  <hwborchers@googlemail.com>
}
\examples{
  x <- 2*runif(5) - 1
  fnTrefethen(x)
  fnWagon(x)

  \dontrun{
  T <- matrix(NA, nrow=1001, ncol=1001)
  for (i in 1:1001) {
    for (j in 1:1001) {
        T[i, j] <- fnTrefethen(c(x[i], y[j]))
    }
  }
  image(x, y, T)
  contour(x, y, T, add=TRUE)
  }
}

\keyword{ testfunctions }