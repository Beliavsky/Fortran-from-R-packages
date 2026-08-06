\name{is.toeplitz}
\alias{is.toeplitz}
\title{ test if argument is a symmetric Toeplitz matrix }
\description{
Auxilary function, used to validate the input of TrenchInverse
}
\usage{
is.toeplitz(x)
}
\arguments{
  \item{x}{ value to be tested }
}
\details{
A symmetric Toeplitz matrix of order n has (i,j)-entry of the
form g[abs(1+i-j)], where g is a vector of length n.
}
\value{ returns True or False according to whether x is or is not
a symmetric Toeplitz matrix
}
\author{ A.I. McLeod }

\seealso{ \code{\link{TrenchInverse}}, \code{\link{toeplitz}} }
\examples{
is.toeplitz(toeplitz(1:5))
is.toeplitz(5)
}
\keyword{ ts }
\keyword{array }

