\name{varFisher}
\alias{varFisher}
\title{
Observed Fisher information matrix
}
\description{
Observed Fisher information matrix on a sample from the type 1 discrete Weibull distribution
}
\usage{
varFisher(x, zero = FALSE)
}
\arguments{
  \item{x}{
a vector of sample values
}
  \item{zero}{
\code{TRUE}, if the support contains \eqn{0}; \code{FALSE} otherwise
}
}
\value{
a list of two matrices: the observed Fisher information matrix, and its inverse, which contains asymptotic variances and covariances of the maximum likelihood estimators of \eqn{q} and \eqn{\beta}}
\author{
Alessandro Barbiero
}
\seealso{\code{\link{estdweibull}}
}
\examples{
x <- rdweibull(100, 2/3, 3/2)
estdweibull(x, "ML")
IF <- varFisher(x)[[2]]
diag(IF)
diag(IF)/length(x)
# asymptotic variances of the ML estimators directly estimated from the sample
}
\keyword{distribution}
\keyword{htest}

