\name{pair.sum}
\alias{pair.sum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pairwise Sums}
\description{
  Computes pairwise sums.
}
\usage{
pair.sum(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric matrix.}
}
\details{
  The function computes all sums of row i and row j with i < j. The function is a wrapper to a C function to do the computation quickly and does no checks concerning the input. 
}
\value{
  Matrix containing the sums.
}

\author{Seija Sirkia}

\seealso{\code{\link{pair.diff}}, \code{\link{pair.prod}}}
\examples{
X <- matrix(1:10, ncol = 2, byrow = FALSE)
pair.sum(X)
}
\keyword{manip }

