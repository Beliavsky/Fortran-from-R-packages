\name{rcov}
\alias{rcov}
\title{Computes robust estimate of covariance matrix}
\usage{
rcov(x)
}
\arguments{
\item{x}{Matrix that you wish to find robust covariance of. Number of
	variables is number of rows, number of observations is number
	of columns. This is the opposite way round to the convention
	expected by \code{var}, for example}
}
\description{
Computes a robust correlation matrix from x.
}
\details{
Method originates from Huber's "Robust Statistics" book.
Note that the columns of \code{x} must be observations, this is the opposite
way around to the usual way for functions like \code{var}.
}
\value{
The robust covariance matrix
}
\seealso{\code{\link{threshold.mwd}}}
\examples{
#
# A standard normal data matrix with 3 variables, 100 observations
#
v <- matrix(rnorm(100*3), nrow=3, ncol=100)
#
# Robust covariance
#
rcov(v)
}
\author{Tim Downie}
\keyword{robust}
\keyword{multivariate}
