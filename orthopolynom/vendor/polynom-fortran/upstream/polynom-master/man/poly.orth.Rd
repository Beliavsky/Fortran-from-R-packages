\name{poly.orth}
\alias{poly.orth}
\title{Construct Orthogonal Polynomials}
\description{
  Construct the orthogonal polynomials on a given vector, up to a
  specified degree.
}
\usage{
poly.orth(x, degree = length(unique(x)) - 1, norm = TRUE)
}
\arguments{
  \item{x}{a numeric vector of abscissae.  When evaluated at \code{x}
    the polynomials will generate an orthonormal set.}
  \item{degree}{maximum degree required.  The default is one fewer than
    the number of distinct values in \code{x}, which is maximum
    possible.}
  \item{norm}{a logical indicating whether the polynomials should be
    normalized.}
}
\value{
  A list of class \code{"polylist"} of objects of class
  \code{"polynomial"} of degree 1, 2, \dots, \code{degree}.
}
\examples{
x <- rep(1:4, 1:4)		# x with repetitions for weighting
x
## [1] 1 2 2 3 3 3 4 4 4 4
polx <- poly.orth(x, 3)		# calculate orthogonal polynomials
polx
## List of polynomials:
## [[1]]
## 0.3162278 
## 
## [[2]]
## -0.9486833 + 0.3162278*x 
## 
## [[3]]
## 2.139203 - 1.863177*x + 0.3450328*x^2 
## 
## [[4]]
## -5.831564 + 8.80369*x - 3.803194*x^2 + 0.4930066*x^3 
v <- sapply(polx, predict, x)   # orthonormal basis
round(crossprod(v), 10)	        # check orthonormality
}
\keyword{symbolmath}
