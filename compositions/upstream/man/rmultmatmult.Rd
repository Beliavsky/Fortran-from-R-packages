\name{rmultmatmult}
\alias{\%*\%.rmult}
\title{inner product for datasets with vector scale}
\description{
An rmult object is considered as a sequence of vectors. The \code{\%*\%}
is considered as the inner multiplication. An inner multiplication with
another vector is the scalar product. an inner multiplication with
a matrix is a matrix multiplication, where the rmult-vectors are either
considered as row or as column vector.  
}
\usage{
\method{\%*\%}{rmult}(x,y)
}
\arguments{
  \item{x}{an rmult vector or dataset of vectors, a numeric vector of
    length (\code{gsi.getD(y)}), or a matrix}
  \item{y}{an rmult vector or dataset of vectors , a numeric vector of
    length (\code{gsi.getD(x)}), or a matrix}
}
\value{
  an object of class \code{"rmult"} or a numeric vector containing the
result of the
corresponding inner products.
}

\details{ The operators try to minic the behavior of \code{\%*\%} on
\code{c()}-vectors as inner product applied in parallel to all vectors of
the dataset. Thus the product of a vector with another \code{rmult}
object or unclassed vector \var{v} results in the scalar product. For
the multiplication with a matrix each vector is considered as a row or
column, whatever is more appropriate. 
}
\note{The product \code{x \%*\% A \%*\% y} is associative. }
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{rmult}}, \code{\link{\%*\%.rmult}}
}

\examples{
x <- rmult(matrix( sqrt(1:12), ncol= 3 ))
x\%*\%x
A <- matrix( 1:9,nrow=3)
x \%*\% A \%*\% x
x \%*\% A
A \%*\% x
x \%*\% 1:3
x \%*\% 1:3
1:3 \%*\% x 

}
\keyword{multivariate}
