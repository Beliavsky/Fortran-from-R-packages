\name{polynomial.functions}
\alias{polynomial.functions}
\title{ Coerce polynomials to functions }
\description{
  This function returns a list with \eqn{n + 1} elements containing the functions of the 
  order $k$ polynomials for orders \eqn{k = 0,\;1,\; \ldots ,\;n}
  and for the given argument \eqn{x}.
}
\usage{
polynomial.functions(polynomials, ...)
}
\arguments{
  \item{polynomials}{a list of polynomial objects}
  \item{\dots}{ further arguments to be passed to or from methods }
}
\details{
  The function uses the method as.function.polynomial to coerce each polynomial object
  to a function object.
}
\value{
  A list of \eqn{n + 1} polynomial objects where each element is the function for the polynomial.
  \item{1 }{Function for the order 0 polynomial}
  \item{2 }{Function for the order 1 polynomial}
  ...
  \item{n+1}{Function for the order \eqn{n} polynomial}
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\examples{
###
### generate a list of T Chebyshev polynomials of
### orders 0 to 10
###
p.list <- chebyshev.t.polynomials( 10, normalized=FALSE )
###
### create the list of functions for each polynomial
###
f.list <- polynomial.functions( p.list )
}
\keyword{ math }
