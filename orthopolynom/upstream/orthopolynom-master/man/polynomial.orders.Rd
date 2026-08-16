\name{polynomial.orders}
\alias{polynomial.orders}
\title{ Create vector of polynomial orders }
\description{
  This function returns a vector with \eqn{n} elements containing the orders of the polynomials
}
\usage{
polynomial.orders(polynomials)
}
\arguments{
  \item{polynomials}{ list of $n$ polynomial objects }
}
\value{
  A vector of \eqn{n} values
  \item{1 }{Order of polynomials[[1]]}
  \item{2 }{Order of polynomials[[2]]}
  ...
  \item{n}{Order of polynomials[[n]]}
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\examples{
###
### generate a list of normalized T Chebyshev polynomials
### of orders 0 to 10
###
p.list <- chebyshev.t.polynomials( 10, normalized=TRUE )
###
### get the vector of polynomial orders
###
p.order <- polynomial.orders( p.list )
}
\keyword{ math }
