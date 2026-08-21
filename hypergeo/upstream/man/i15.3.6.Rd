\name{i15.3.6}
\alias{i15.3.6}
\alias{i15.3.7}
\alias{i15.3.8}
\alias{i15.3.9}
\alias{j15.3.6}
\alias{j15.3.7}
\alias{j15.3.8}
\alias{j15.3.9}
\title{Helper functions}
\description{
Helper functions for equations 15.3.6-15.3.9
}
\usage{
i15.3.6(A, B, C)
i15.3.7(A, B, C)
i15.3.8(A, B, C)
i15.3.9(A, B, C)
j15.3.6(A, B, C)
j15.3.7(A, B, C)
j15.3.8(A, B, C)
j15.3.9(A, B, C)
}
\arguments{
  \item{A,B,C}{Parameters of the hypergeometric function}
}
\details{

These functions are named for their equation numbers in Abramowitz and
Stegun.
  
Functions \code{i15.3.?()} return the factors at the beginning of
equations 15.3.6-9.  These functions return zero if the denominator is
infinite (because it includes a gamma function of a nonpositive
integer). 


Functions \code{j15.3.?()} check for the appropriate
arguments of the gamma function being nonpositive integers.

}
\references{
  M. Abramowitz and I. A. Stegun 1965. \emph{Handbook of
    mathematical functions}.  New York: Dover
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{hypergeo}}}
\examples{
i15.3.6(1.1, 3.2, pi)
}
\keyword{math}
