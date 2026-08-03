\name{kurtosis}
\alias{kurtosis}
\title{ Pearson's measure of kurtosis }
\description{
  This function computes the estimator of Pearson's measure of kurtosis.
}
\usage{
kurtosis(x, na.rm = FALSE)
}
\arguments{
  \item{x}{ a numeric vector, matrix or data frame. }
  \item{na.rm}{ logical. Should missing values be removed? }
}
\author{ Lukasz Komsta }
\seealso{ \code{\link{geary}}, \code{\link{anscombe.test}} }
\examples{
set.seed(1234)
kurtosis(rnorm(1000))
}
\keyword{univar}

