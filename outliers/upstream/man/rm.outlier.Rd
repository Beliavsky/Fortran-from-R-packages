\name{rm.outlier}
\alias{rm.outlier}
\title{ Remove the value(s) most differing from the mean }
\description{

If the outlier is detected and confirmed by statistical tests, this function
can remove it or replace by sample mean or median.

}
\usage{
rm.outlier(x, fill = FALSE, median = FALSE, opposite = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a dataset, most frequently a vector. If argument is a dataframe, then outlier is removed 
from each column by \code{sapply}. The same behavior is applied by \code{apply} when the matrix 
is given.
 }
  \item{fill}{ If set to TRUE, the median or mean is placed instead of outlier. Otherwise, 
	the outlier(s) is/are simply removed. }
  \item{median}{ If set to TRUE, median is used instead of mean in outlier replacement. }
  \item{opposite}{ if set to TRUE, gives opposite value (if largest value has maximum difference
from the mean, it gives smallest and vice versa)}
}

\value{

A dataset of the same type as argument, with outlier(s) removed or replacement by appropriate means
or medians.

}
\author{ Lukasz Komsta}
\seealso{ \code{\link{outlier}} }
\examples{
set.seed(1234)
y=rnorm(100)
outlier(y)
outlier(y,opposite=TRUE)
rm.outlier(y)
rm.outlier(y,opposite=TRUE)
dim(y) <- c(20,5)
outlier(y)
outlier(y,logical=TRUE)
outlier(y,logical=TRUE,opposite=TRUE)
rm.outlier(y)
rm.outlier(y,opposite=TRUE)
}
\keyword{htest}
