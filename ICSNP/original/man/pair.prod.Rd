\name{pair.prod}
\alias{pair.prod}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Pairwise Products }
\description{
  Computes pairwise elementwise products.
}
\usage{
pair.prod(X)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{a numeric matrix.}
}
\details{
  The function computes all elementwise products of row i and row j with i < j. The function is a wrapper to a C function to do the computation quickly and does no checks concerning the input. 
}
\value{
  Matrix containing the products.
}
\author{Klaus Nordhausen}
\seealso{\code{\link{pair.diff}}, \code{\link{pair.sum}}}
\examples{
X <- matrix(1:10, ncol = 2, byrow = FALSE)
pair.prod(X)
}
\keyword{manip }
