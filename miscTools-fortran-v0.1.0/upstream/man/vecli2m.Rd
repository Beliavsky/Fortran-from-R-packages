\name{vecli2m}
\alias{vecli2m}
\title{Convert vector of linear independent values into a Matrix}

\description{
   Converts a vector into a symmetric matrix that the original vector
   contains the linear independent values of the returned symmetric
   matrix.
}

\usage{vecli2m( v )}

\arguments{
   \item{ v }{a vector.}
}

\seealso{\code{\link{vecli}}, \code{\link{veclipos}}.}

\author{Arne Henningsen}

\examples{
   v <- c( 11, 12, 13, 22, 23, 33 )
   vecli2m( v )
}

\keyword{array}
