\name{matmult}
\alias{\%*\%}
\alias{\%*\%.default}
\title{inner product for matrices and vectors}
\description{
  Multiplies two matrices, if they are conformable. If one argument
     is a vector, it will be coerced to either a row or a column matrix
     to make the two arguments conformable. If both are vectors it will
     return the inner product.  
   }
\usage{
x \%*\% y
\method{\%*\%}{default}(x, y)
}
\arguments{
  \item{x,y}{numeric or complex matrices or vectors}
}
\value{
  The matrix product.  Uses 'drop' to get rid of dimensions which
  have only one level.
}

\details{ 
  This is a copy of the %\code{\link[base]{\%*\%}} 
  \code{base::\%*\%}
  function. The
  function is made generic to allow the definition of specific methods.
}
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{\%*\%.rmult}}
}
\examples{
M <- matrix(c(
0.2,0.1,0.0,
0.1,0.2,0.0,
0.0,0.0,0.2),byrow=TRUE,nrow=3)
x <- c(1,1,2)
M \%*\% x
x \%*\% M
x \%*\% x
M \%*\% M
t(x) \%*\% M

}
\keyword{multivariate}
