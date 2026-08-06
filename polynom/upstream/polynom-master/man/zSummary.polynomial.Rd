\name{zSummary.polynomial}
\alias{Summary.polynomial}
\title{Summary Group Methods for Polynomials}
\description{
  Allows summary group generics to be used on polynomial arguments.
}
\usage{
\method{Summary}{polynomial}(\dots, na.rm = FALSE)
}
\arguments{
  \item{\dots}{\R objects, the first supplied of class
    \code{"polynomial"}.} 
  \item{na.rm}{logical: should missing values be removed?}
}
\details{
  For the \code{sum} and \code{prod} functions, the sum and product of
  the given polynomials, respectively.  For the other members of the
  Summary group, an error is returned.

  Note that one could \emph{order} polynomials by divisibility, and
  define \code{min} and \code{max} as the corresponding lattice meet and
  join, i.e., the greatest common divisor and the least common multiple,
  respectively.  This is currently not provided: instead, functions
  \code{\link{GCD}} and \code{\link{LCM}} should be called directly.
}
\seealso{
  \code{\link{Math.polynomial}},
  \code{\link{Ops.polynomial}}
}
\keyword{symbolmath}
