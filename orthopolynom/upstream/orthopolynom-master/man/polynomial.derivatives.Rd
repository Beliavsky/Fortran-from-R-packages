\name{polynomial.derivatives}
\alias{polynomial.derivatives}
\title{ Create list of polynomial derivatives }
\description{
  This function returns a list with \eqn{n + 1} elements containing polynomial objects 
  which are the derivatives of the order \eqn{k} polynomials 
  for orders \eqn{k = 0,\;1,\; \ldots ,\;n}.
}
\usage{
polynomial.derivatives(polynomials)
}
\arguments{
  \item{polynomials}{ list of polynomial objects }
}
\details{
  The polynomial objects in the argument \code{polynomials} are as follows
  \itemize{
  \item{1}{order 0 polynomial}
  \item{2}{order 1 polynomial}
  ...
  \item{n+1}{order \eqn{n} polynomial}
  }
}
\value{
  List of \eqn{n + 1} polynomial objects
  \item{1 }{derivative of \code{polynomials[[1]]}}
  \item{2 }{derivative of \code{polynomials[[2]]}}
  ...
  \item{n+1 }{derivative of \code{polynomials[[n+1]]}}
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\examples{
###
### generate a list of normalized T Chebyshev polynomials of
### orders 0 to 10
###
p.list <- chebyshev.t.polynomials( 10, normalized=TRUE )
###
### generate the corresponding list of polynomial derivatives
###
p.deriv <- polynomial.derivatives( p.list )
}
\keyword{ math }
