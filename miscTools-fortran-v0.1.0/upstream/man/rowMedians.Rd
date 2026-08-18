\name{rowMedians}
\alias{rowMedians}
\title{Medians of Rows}

\description{
   Compute the sample medians of the rows of a data.frame or matrix.
}

\usage{
rowMedians( x, na.rm = FALSE )
}

\arguments{
   \item{x}{a data.frame or matrix.}
   \item{na.rm}{a logical value indicating whether \code{NA} values should be
      stripped before the computation proceeds.}
}

\value{
   A vector of the medians of each row of \code{x}.
}

\seealso{\code{\link{colMedians}},\code{\link{median}},\code{\link{colMeans}}.}

\author{Arne Henningsen}

\examples{
   m <- matrix( 1:12, nrow = 4 )
   rowMedians( m )
}

\keyword{array}
