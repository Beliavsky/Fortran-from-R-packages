\name{rowSums.dist}
\alias{rowSums.dist}
\alias{colSums.dist}
\alias{rowMeans.dist}
\alias{colMeans.dist}
\title{Row Sums/Means of Sparse Symmetric Matrices}
\description{
  Compute the row (column) sums or means for a sparse symmetric 
  (distance) matrix.
}
\usage{
rowSums.dist(x, na.rm = FALSE)
rowMeans.dist(x, na.rm = FALSE, diag = TRUE)

colSums.dist(x, na.rm = FALSE)
colMeans.dist(x, na.rm = FALSE, diag = TRUE)
}
\arguments{
  \item{x}{an object of class \code{dist}.}
  \item{na.rm}{logical, should missing values (including \code{NaN}) be
    omitted from the summation?}
  \item{diag}{logical, should the diagonal elements be included
    in the computation?}
}
\details{
  These functions are more efficient than expanding an object of
  class \code{dist} to matrix and using \code{rowSums} or \code{rowMeans}.

  \code{colSums} and \code{colMeans} are provided for convenience.
  However, note that due to symmetry the result is always the
  same as for \code{rowSums} or \code{rowMeans}.
}
\value{
  A numeric vector of row sums.
}
\author{Christian Buchta}
\seealso{\code{as.matrix}, \code{as.dist}, and \code{rowSums}.}
\examples{
## 
x <- matrix(runif(10*2),ncol=2)
d <- dist(x)
rowSums(as.matrix(d))
rowSums.dist(d)			# the same

rowMeans(as.matrix(d))
rowMeans.dist(d)		# the same
rowMeans.dist(d, diag = FALSE)	# not the same
## NAs
d[3] <- NA
rowSums.dist(d, na.rm = TRUE)
rowMeans.dist(d, na.rm = TRUE)
}
\keyword{cluster}
