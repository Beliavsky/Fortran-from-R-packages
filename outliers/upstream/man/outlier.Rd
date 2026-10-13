\name{outlier}
\alias{outlier}
\title{ Find value with largest difference from the mean }
\description{
  Finds value with largest difference between it and sample mean, which can be an outlier.
}
\usage{
outlier(x, opposite = FALSE, logical = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a data sample, vector in most cases. If argument is a dataframe, then outlier is calculated for
each column by \code{sapply}. The same behavior is applied by \code{apply} when the matrix is given.
  }
  \item{opposite}{ if set to TRUE, gives opposite value (if largest value has maximum difference
from the mean, it gives smallest and vice versa)}
  \item{logical}{ if set to TRUE, gives vector of logical values, and possible outlier position
 is marked by TRUE}
}
\value{

A vector of value(s) with largest difference from the mean.

}
\author{ Lukasz Komsta, corrections by Markus Graube }
\seealso{ \code{\link{rm.outlier}}}
\examples{
set.seed(1234)
y=rnorm(100)
outlier(y)
outlier(y,opposite=TRUE)
dim(y) <- c(20,5)
outlier(y)
outlier(y,opposite=TRUE)
}
\keyword{htest}
