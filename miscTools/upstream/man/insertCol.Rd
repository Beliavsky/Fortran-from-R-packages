\name{insertCol}
\alias{insertCol}
\title{Insert Column into a Matrix}

\description{
   Insert a new column into a matrix.
}

\usage{insertCol( m, c, v = NA, cName = "" )}

\arguments{
   \item{ m }{matrix.}
   \item{ c }{column number where the new column should be inserted.}
   \item{ v }{optional values of the new column.}
   \item{ cName }{optional character string: the name of the new column.}
}

\value{
   a matrix with one more column than the provided matrix \code{m}.
}

\seealso{\code{\link{insertRow}}.}

\author{Arne Henningsen}

\examples{
   m <- matrix( 1:4, 2 )
   insertCol( m, 2, 5:6 )
}

\keyword{array}
