\name{gsi2.invperm}
\alias{gsi2.invperm}
\title{Internal function: Invert a permutation}
\description{
  Finds the inverse of a permutation given as a vector of indices.
}
\usage{
gsi2.invperm( i,n )
}
\arguments{
  \item{i}{ a sequence of different integers in \code{1:n} considered as
  a permutation given by \code{p=unique(c(i,1:n))}}
  \item{n}{ the number of elements to be permuted}
}
\value{
  an integer vector \var{v} describing the inverse permutation of p.
}
\note{Do not use gsi.* functions directly since they are internal
  functions of the package}
\details{
  The inverse permutation is defined by \code{p[\var{v}]==1:n} and
  \code{v[\var{p}]==1:n}.  
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
    \code{\link{gsi}}
}
\examples{
#gsi2.invperm(c(2,3),10)
}
\keyword{internal}
