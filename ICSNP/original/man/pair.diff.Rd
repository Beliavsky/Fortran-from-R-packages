\name{pair.diff}
\alias{pair.diff}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Pairwise Differences}
\description{
  Computes pairwise differences.
}
\usage{
pair.diff(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric matrix.}
}
\details{
  The function computes all differences of row i and row j with i < j. The function is a wrapper to a C function to do the computation quickly and does no checks concerning the input.
}
\value{
  Matrix containing the differences.
}
\author{Seija Sirkia}
\seealso{\code{\link{pair.prod}}, \code{\link{pair.sum}}}
\examples{
X <- matrix(1:10, ncol = 2, byrow = FALSE)
pair.diff(X)
}
\keyword{manip }
