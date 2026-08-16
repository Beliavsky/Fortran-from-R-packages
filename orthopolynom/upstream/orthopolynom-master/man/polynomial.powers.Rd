\name{polynomial.powers}
\alias{polynomial.powers}
\title{ Create a list of polynomial linear combinations }
\description{
  This function returns a list with \eqn{n + 1} elements containing the vector of 
  linear combinations of the order \eqn{k} polynomials for orders \eqn{k = 0,\;1,\; \ldots ,\;n}.
  Each element in the list is a vector.
}
\usage{
polynomial.powers(polynomials)
}
\arguments{
  \item{polynomials}{ A list of polynomials }
}
\details{
  The \eqn{j}-th component in the list is a vector of linear combinations of
  the order \eqn{k} polynomials for orders \eqn{k = 0,\;1,\; \ldots ,\;j - 1} equal to
  the monomial x raised to the power \eqn{j - 1}.
}
\value{
  A list of \eqn{n + 1} elements where each element is a vector of linear combinations.
  \item{1 }{Linear combination(s) of polynomials up to order 0}
  \item{2 }{Linear combination(s) of polynomials up to order 1}
  ...
  \item{n+1}{Linear combination(s) of polynomials up to order \eqn{n}}
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\examples{
###
### generate Legendre polynomials of orders 0 to 10
###
polynomials <- legendre.polynomials( 10 )
###
### generate list of linear combinations of these polynomials
###
alphas <- polynomial.powers( polynomials )
print( alphas )
}
\keyword{ math }
