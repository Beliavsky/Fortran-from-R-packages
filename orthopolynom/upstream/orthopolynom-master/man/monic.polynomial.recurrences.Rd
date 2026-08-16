\name{monic.polynomial.recurrences}
\alias{monic.polynomial.recurrences}
\title{ Create data frame of monic recurrences }
\description{
  This function returns a data frame with parameters required to construct monic orthogonal polynomials 
  based on the standard recurrence relation for the non-monic polynomials.  
  The recurrence relation for monic orthogonal polynomials is as follows.
  \deqn{q_{k + 1} \left( x \right) = \left( {x - a_k } \right)\;q_k \left( x \right) - b_k \;q_{k - 1} \left( x \right)}
  We require that \eqn{q_{-1} \left( x \right) = 0} and \eqn{q_0 \left( x \right) = 1}.  
  The recurrence for non-monic orthogonal polynomials is given by
  \deqn{c_k \;p_{k + 1} \left( x \right) = \left( {d_k  + e_k \;x} \right)\;p_k \left( x \right) - f_k \;p_{k - 1} \left( x \right)}
  We require that \eqn{p_{-1} \left( x \right) = 0} and \eqn{p_0 \left( x \right) = 1}.  
  The monic polynomial recurrence parameters, \bold{a} and \bold{b}, are related to 
  the non-monic polynomial parameter vectors \bold{c}, \bold{d}, \bold{e} and \bold{f} in the following manner.
  \deqn{a_k  =  - \frac{{d_k }}{{e_k }}}
  \deqn{b_k  = \frac{{c_{k - 1} \;f_k }}{{e_{k - 1} \;e_k }}} with \eqn{b_0 = 0}.
}
\usage{
monic.polynomial.recurrences(recurrences)
}
\arguments{
  \item{recurrences}{ the data frame of recurrence parameter vectors \bold{c}, \bold{d}, \bold{e} and \bold{f} }
}
\value{
  A data frame with \eqn{n + 1} rows and two named columns, \bold{a} and \bold{b}.
}
\references{
Abramowitz, M. and I. A. Stegun, 1968. \emph{Handbook of Mathematical Functions with
Formulas, Graphs, and Mathematical Tables}, Dover Publications, Inc., New York.

Courant, R., and D. Hilbert, 1989. \emph{Methods of Mathematical Physics},
John Wiley, New York, NY.

Press, W. H., S. A. Teukolsky, W. T. Vetterling, and B. P. Flannery, 1992. 
\emph{Numerical Recipes in C}, Cambridge University Press, Cambridge, U.K.

Szego, G., 1939. \emph{Orthogonal Polynomials}, 23, American Mathematical Society
Colloquium Publications, Providence, RI.
}
\author{ Frederick Novomestky \email{fnovomes@poly.edu} }
\seealso{ 
  \code{\link{orthogonal.polynomials}},
}
\examples{
###
### construct a list of the recurrences for the T Chebyshev polynomials of
### orders 0 to 10
###
r <- chebyshev.t.recurrences( 10, normalized=TRUE )
###
### construct the monic polynomial recurrences from the above list
###
m.r <- monic.polynomial.recurrences( r )
}
\keyword{ math }
