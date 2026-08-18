\name{insertRow}
\alias{insertRow}
\title{Insert Row into a Matrix}

\description{
   Insert a new row into a matrix.
}

\usage{insertRow( m, r, v = NA, rName = "" )}

\arguments{
   \item{ m }{matrix.}
   \item{ r }{row number where the new row should be inserted.}
   \item{ v }{optional values for the new row.}
   \item{ rName }{optional character string: the name of the new row.}
}

\value{
   a matrix with one more row than the provided matrix \code{m}.
}

\seealso{\code{\link{insertCol}}.}

\author{Arne Henningsen}

\examples{
   m <- matrix( 1:4, 2 )
   insertRow( m, 2, 5:6 )
}

\keyword{array}
