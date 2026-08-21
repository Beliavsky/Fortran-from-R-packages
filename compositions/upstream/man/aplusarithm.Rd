\name{aplusarithm}
\alias{+.aplus}
\alias{-.aplus}
\alias{*.aplus}
\alias{/.aplus}
\alias{perturbe.aplus}
\alias{power.aplus}
\title{vectorial arithmetic for data sets with aplus class}
\description{
  The positive vectors equipped with the perturbation (defined as
  the element-wise product) as Abelian sum, and powertransform (defined as the element-wise
 powering with a scalar) as scalar multiplication forms a real vector
 space. These vector space operations are defined here in a similar way
 to \code{\link{+.rmult}}.
}
\usage{
perturbe.aplus(x,y)
\method{+}{aplus}(x,y)
\method{-}{aplus}(x,y)
\method{*}{aplus}(x,y)
\method{/}{aplus}(x,y)
##  Methods for aplus
##   x+y
##   x-y
##   -x
##   x*r
##   r*x
##   x/r
power.aplus(x,r)
}
\arguments{
  \item{x}{an aplus vector or dataset of vectors }
  \item{y}{an aplus vector or dataset of vectors }
  \item{r}{a numeric vector of size 1 or nrow(x)}
}
\value{
  an object of class \code{"aplus"} containing the result of the
  corresponding operation on the vectors.  
}
\details{
  
  
  The operators try to mimic the parallel operation of R for vectors of
real numbers to vectors of amounts, represented as matrices containing
the vectors as rows and works like the operators for \code{{\link{rmult}}} 

}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{rmult}}, \code{\link{\%*\%.rmult}}
}

\examples{
x <- aplus(matrix( sqrt(1:12), ncol= 3 ))
x
x+x
x + aplus(1:3)
x * 1:4
1:4 * x
x / 1:4
x / 10
power.aplus(x,1:4)
}
\keyword{multivariate}
