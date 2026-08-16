\name{monic.polynomials}
\alias{monic.polynomials}
\title{ Create list of monic orthogonal polynomials }
\description{
  This function returns a list with \eqn{n + 1} elements containing 
  the order \eqn{k} monic polynomials 
  for orders \eqn{k = 0,\;1,\; \ldots ,\;n}.
}
\usage{
monic.polynomials(monic.recurrences)
}
\arguments{
  \item{monic.recurrences}{ a data frame containing the parameters a and b }
}
\value{
  A list with \eqn{n + 1} polynomial objects
  \item{1 }{order 0 monic orthogonal polynomial}
  \item{2 }{order 1 monic orthogonal polynomial}
  ...
  \item{n+1 }{order \eqn{n} monic orthogonal polynomial}
}
\references{
Press, W. H., S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery, 1992. 
\emph{Numerical Recipes in C}, Cambridge University Press, Cambridge, U.K.
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{
  \code{\link{monic.polynomial.recurrences}}
}
\examples{
###
### generate the recurrences for the T Chebyshev polynomials
### of orders 0 to 10
###
r <- chebyshev.t.recurrences( 10, normalized=TRUE )
###
### get the corresponding monic polynomial recurrences
###
m.r <- monic.polynomial.recurrences( r )
###
### obtain the list of monic polynomials
###
p.list <- monic.polynomials( m.r )
}
\keyword{ math }
