\name{veclipos}
\alias{veclipos}
\title{Position in a vector of linear independent values}

\description{
   Returns the position of the [\code{i},\code{j}]th element of a symmetric
   \code{n} \eqn{\times}{x} \code{n} matrix that this element has in a vector
   of the linear independent values of the matrix.
}

\usage{veclipos( i, j, n )}

\arguments{
   \item{i}{row of the element in the matrix.}
   \item{j}{column of the element in the matrix.}
   \item{n}{dimension of the matrix.}
}

\note{
   A symmetric \code{n} \eqn{\times}{x} \code{n} matrix has n*(n+1)/2
   independent values.\cr
   The function is: n*(n-1)/2-((n-min(i,j))*(n-min(i,j)+1)/2)+max(i,j)
}

\seealso{\code{\link{vecli}}, \code{\link{vecli2m}}.}

\author{Arne Henningsen}

\examples{
   veclipos( 1, 2, 3 )  # returns: 2
}

\keyword{array}
