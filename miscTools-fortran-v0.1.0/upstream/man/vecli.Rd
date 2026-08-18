\name{vecli}
\alias{vecli}
\title{Vector of linear independent values}

\description{
   Returns a vector containing the linear independent elements of a symmetric
   matrix (of full rank).
}

\usage{vecli( m )}

\arguments{
   \item{ m }{symmetric matrix}
}

\seealso{\code{\link{veclipos}}.}

\author{Arne Henningsen}

\examples{
   # a symmetric n x n matrix
   m <- cbind(c(11,12,13),c(12,22,23),c(13,23,33))
   vecli(m)  # returns: 11 12 13 22 23 33
}

\keyword{array}
