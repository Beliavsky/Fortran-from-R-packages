\name{jacobi.matrices}
\alias{jacobi.matrices}
\title{ Create list of Jacobi matrices from monic recurrence parameters }
\description{
  Return a list of $n$ real symmetric, tri-diagonal matrices which are 
  the principal minors of the \eqn{n \times n} Jacobi matrix derived from 
  the monic recurrence parameters, \eqn{a} and \eqn{b}, for orthogonal polynomials.
}
\usage{
jacobi.matrices(r)
}
\arguments{
  \item{r}{ a data frame containing the parameters \eqn{a} and \eqn{b} }
}
\value{
  A list of symmetric, tri-diagnonal matrices
  \item{1 }{a \eqn{1 \times 1} matrix}
  \item{2 }{a \eqn{2 \times 2} matrix}
  ...
  \item{n }{an \eqn{n \times n} matrix}
}
\references{
Press, W. H., S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery, 1992. 
\emph{Numerical Recipes in C}, Cambridge University Press, Cambridge, U.K.
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\examples{
r <- chebyshev.t.recurrences( 5 )
m.r <- monic.polynomial.recurrences( r )
j.m <- jacobi.matrices( m.r )
}
\keyword{ math }
