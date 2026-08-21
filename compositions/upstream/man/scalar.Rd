\name{scalar}
\alias{scalar}
\alias{scalar.default}
\title{Parallel scalar products }
\description{
Computes scalar products of datasets of vectors or vectorial quantities.
}
\usage{
scalar(x,y)
\method{scalar}{default}(x,y)
          }
\arguments{
  \item{x}{a vector or a matrix with rows considered as vectors}
  \item{y}{a vector or a matrix with rows considered as vectors}
}
\value{
  a numerical vector containing the scalar products of the vectors given
  by x and y. If both \code{x} and \code{y} contain more than one
  vector the function uses parallel operation like it would happen with
  an ordinary product of vectors.
}
\details{
  The scalar product of two vectors is defined as:
  \deqn{scalar(x,y):=  \sum(x_iy_i) }
}
\note{
  The computation of the scalar product implicitly applies
  the \code{\link{cdt}} transform, which implies that the scalar products
  corresponding to the given geometries are returned for \code{\link{acomp}},
  \code{\link{rcomp}}, \code{\link{aplus}},
  \code{\link{rplus}}-objects. Even a useful scalar product for factors
  is induced in this way.

}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\examples{
scalar(acomp(c(1,2,3)),acomp(c(1,2,3)))
scalar(rmult(c(1,2,3)),rmult(c(1,2,3)))
}
\keyword{multivariate}
