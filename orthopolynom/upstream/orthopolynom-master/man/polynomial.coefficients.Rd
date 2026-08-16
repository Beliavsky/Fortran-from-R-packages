\name{polynomial.coefficients}
\alias{polynomial.coefficients}
\title{ Create list of polynomial coefficient vectors }
\description{
  This function returns a list with \eqn{n + 1} elements containing 
  the vector of coefficients of the order \eqn{k} polynomials 
  for orders \eqn{k = 0,\;1,\; \ldots ,\;n}.
  Each element in the list is a vector.
}
\usage{
polynomial.coefficients(polynomials)
}
\arguments{
  \item{polynomials}{ list of polynomial objects }
}
\value{
  A list of \eqn{n + 1} polynomial objects where each element is a vector of coefficients.
  \item{1 }{Coefficient(s) of the order 0 polynomial}
  \item{2 }{Coefficient(s) of the order 1 polynomial}
  ...
  \item{n+1}{Coefficient(s) of the order \eqn{n} polynomial}
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\examples{
###
### generate a list of normalized T Chebyshev polynomials
### of orders 0 to 10
###
p.list <- chebyshev.t.polynomials( 10, normalized=TRUE )
###
### obtain the list of coefficients for these polynomials
###
p.coef <- polynomial.coefficients( p.list )
}
\keyword{ math }
