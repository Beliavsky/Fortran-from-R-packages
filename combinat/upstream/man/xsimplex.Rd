\name{xsimplex}
\alias{xsimplex}
\title{  Generates all points on a (p,n) simplex lattice (i.e. a p-part
             composition of n).}
\description{
             Generates all points on a {p,n} simplex lattice (i.e. a p-part
             composition of n).  Each point is represented as x, a
             p-dimensional vector of nonnegative integers that sum to n.
             If argument "fun" is not null,  applies a function given
             by the argument to each point.  If simplify is FALSE,  returns
             a list; else returns a vector or an array.  "..." are passed
             unchanged to function given by argument fun,  if any. }
\usage{
xsimplex(p, n, fun=NULL, simplify=TRUE, ...)
}
%- maybe also `usage' for other objects documented here.
\arguments{
  \item{p}{ first parameter of lattice description }
  \item{n}{  second parameter of lattice description}
  \item{fun}{ function to be applied pointwise }
  \item{simplify}{ logical: if FALSE, value is a list, otherwise a vector or
 array}
  \item{\dots}{ parameters to be passed to \code{fun} }
}
%\details{
%}
%\value{
%}
%\references{ }
%
\examples{
#Compute Multinomial(n = 4, pi = rep(1/3, 3)) p.f.:
   xsimplex(3, 4, dmnom, prob=1/3)
}
\keyword{models}

