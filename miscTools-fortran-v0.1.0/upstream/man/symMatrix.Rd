\name{symMatrix}
\alias{symMatrix}
\title{Symmetric Matrix}

\description{
   Create a Symmetric Matrix.
}

\usage{
   symMatrix( data = NA, nrow = NULL, byrow = FALSE,
   upper = FALSE )
}

\arguments{
   \item{data}{an optional data vector.}
   \item{nrow}{the desired number of rows and columns.}
   \item{byrow}{logical. If 'FALSE' (the default) the matrix is filled by
      columns, otherwise the matrix is filled by rows.}
   \item{upper}{logical. If 'FALSE' (the default) the lower triangular part
      of the matrix (including the diagonal) is filled, otherwise 
      the upper triangular part of the matrix is filled.}
}


\value{
   a symmetric matrix.
}

\seealso{\code{\link{matrix}}, \code{\link{lower.tri}}.}

\author{Arne Henningsen}

\examples{
   # fill the lower triangular part by columns
   symMatrix( 1:10, 4 )
   # fill the upper triangular part by columns
   symMatrix( 1:10, 4, upper = TRUE )
   # fill the lower triangular part by rows
   symMatrix( 1:10, 4, byrow = FALSE )
}

\keyword{array}
