\name{triang}
\alias{triang}
\title{Upper triangular matrix from a vector}

\description{
   Creates an upper triangular square matrix from a vector.
}

\usage{triang( v, n )}

\arguments{
   \item{v}{vector}
   \item{n}{desired dimension of the returned square matrix}
}

\note{
   If the vector has less elements than the upper triangular matrix,
   the last elements are set to zero.
}

\seealso{\code{\link{veclipos}}.}

\author{Arne Henningsen}

\examples{
   v <- c( 1:5 )
   triang( v, 3 )
}

\keyword{array}
