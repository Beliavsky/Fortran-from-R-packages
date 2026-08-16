\name{polynomial.values}
\alias{polynomial.values}
\title{ Create vector of polynomial values }
\description{
  This function returns a list with \eqn{n + 1} elements containing the values of the 
  order \eqn{k} polynomials for orders \eqn{k = 0,\;1,\; \ldots ,\;n}
  and for the given argument \eqn{x}.
}
\usage{
polynomial.values( polynomials, x )
}
\arguments{
  \item{polynomials}{ list of polynomial objects }
  \item{x}{ the argument which can be any numeric object }
}
\value{
  A list of \eqn{n + 1} polynomial objects where each element is the value of the polynomial.
  \item{1 }{Value(s) for the order 0 polynomial}
  \item{2 }{Value(s) for the order 1 polynomial}
  ...
  \item{n+1}{Value(s) for the order \eqn{n} polynomial}
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\examples{
###
### generate a list of T Chebyshev polynomials of
### orders 0 to 10
###
p.list <- chebyshev.t.polynomials( 10, normalized=FALSE )
x <- seq( -2, 2, .01 )
###
### compute the value of the polynomials for the given range of values in x
###
y <- polynomial.values( p.list, x )
print( y )
}
\keyword{ math }
