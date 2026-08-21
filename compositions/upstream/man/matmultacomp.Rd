\name{acompscalarproduct}
\alias{\%*\%.acomp}
\alias{\%*\%.aplus}
\title{inner product for datasets with a vector space structure}
\description{  
acomp and aplus objects are considered as (sets of) vectors. The
\code{\%*\%} is considered as the inner multiplication. An inner
multiplication with another vector is the scalar product. An inner
multiplication with a matrix is a matrix multiplication, where the
vectors are either considered as row or as column vector.

}
\usage{
\method{\%*\%}{acomp}(x,y)
\method{\%*\%}{aplus}(x,y)
}
\arguments{
  \item{x}{a acomp or aplus object or a matrix interpreted in clr, ilr or ilt coordinates}
  \item{y}{a acomp or aplus object or a matrix interpreted in clr, ilr or ilt coordinates}
}
\value{
  Either a numeric vector containing the scalar products, or an object of
  type acomp or aplus containing the vectors transformed with the given
  matrix. 
}

\details{ The operators try to mimic the behavior of \code{\%*\%} on
\code{c()}-vectors as inner product, applied in parallel to all row-vectors of
the dataset. Thus the product of a vector with a vector of the same
type results in the scalar product of both. For the multiplication with a matrix
each vector is considered as a row or column, whatever is more
appropriate. The matrix itself is considered as representing a linear
mapping (endomorphism) of the vector space to a space of the same type. The mapping is
represented in clr, ilr or ilt coordinates. Which of the aforementioned
coordinate systems is used is judged from the type of \var{x} and from
the dimensions of the \var{A}. 
} 
\author{K.Gerald v.d. Boogaart \url{http://www.stat.boogaart.de}}
\seealso{
  \code{\link{\%*\%.rmult}}
}

\examples{
x <- acomp(matrix( sqrt(1:12), ncol= 3 ))
x\%*\%x
A <- matrix( 1:9,nrow=3)
x \%*\% A \%*\% x
x \%*\% A
A \%*\% x
A <- matrix( 1:4,nrow=2)
x \%*\% A \%*\% x
x \%*\% A
A \%*\% x
x <- aplus(matrix( sqrt(1:12), ncol= 3 ))
x\%*\%x
A <- matrix( 1:9,nrow=3)
x \%*\% A \%*\% x
x \%*\% A
A \%*\% x
}
\keyword{multivariate}
