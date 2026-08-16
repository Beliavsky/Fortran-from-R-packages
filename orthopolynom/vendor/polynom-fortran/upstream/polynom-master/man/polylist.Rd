\name{polylist}
\alias{polylist}
\alias{as.polylist}
\alias{is.polylist}
\alias{[.polylist}
\alias{c.polylist}
\alias{deriv.polylist}
\alias{integral.polylist}
\alias{plot.polylist}
\alias{print.polylist}
\alias{rep.polylist}
\alias{unique.polylist}
\title{Lists of Polynomials}
\description{Create and manipulate lists of polynomials.}
\usage{
polylist(\dots)
as.polylist(x)
is.polylist(x)
}
\arguments{
  \item{\dots}{\R objects, polynomials or capable of coercion to polynomials.}
  \item{x}{an \R object.}
}
\details{
  \code{polylist} takes a list of arguments, tries to convert each into
  a polynomial (see \code{\link{polynomial}}), and sets the class of the
  list to \code{"polylist"}.

  \code{as.polylist} tries to coerce its arguments to a polylist, and
  will do so for arguments which are polynomials or lists thereof.

  \code{is.polylist} tests whether its argument is a polylist.

  This class has several useful methods, such as taking derivatives
  (\code{\link{deriv}}) and antiderivatives (\code{\link{integral}}),
  printing and plotting, subscripting, computing sums and products of
  the elements, and methods for \code{\link{c}}, \code{\link{rep}}, and
  \code{\link{unique}}.
}
\examples{
## Calculate orthogonal polynomials
pl <- poly.orth(rep(1:4, 1:4), 3)
pl
plot(pl)
deriv(pl)
integral(pl)
sum(pl)
prod(pl)
unique(rep(pl, 3)[c(8, 12)])
}
\keyword{symbolmath}
