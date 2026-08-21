\name{rmultarithm}
\alias{+.rmult}
\alias{-.rmult}
\alias{*.rmult}
\alias{/.rmult}
\title{vectorial arithmetic for datasets in a classical vector scale}
\description{
  vector space operations computed for multiple vectors in parallel
}
\usage{
## Methods for class "rmult"
## x+y
## x-y
## -x
## x*r
## r*x
## x/r
}
\arguments{
  \item{x}{an rmult vector or dataset of vectors }
  \item{y}{an rmult vector or dataset of vectors }
  \item{r}{a numeric vector of size 1 or nrow(x)}
}
\value{
  an object of class \code{"rmult"} containing the result of the
  corresponding operation on the vectors.  
}
\details{
The operators try to mimic the parallel operation of R on vectors of
real numbers on vectors of vectors represented as matrices containing
the vectors as rows. 

}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{rmult}}, \code{\link{\%*\%.rmult}}
}

\examples{
x <- rmult(matrix( sqrt(1:12), ncol= 3 ))
x
x+x
x + rmult(1:3)
x * 1:4
1:4 * x
x / 1:4
x / 10
}
\keyword{multivariate}
