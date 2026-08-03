\name{geary}
\alias{geary}
\title{ Geary's measure of kurtosis }
\description{
  This function computes an estimator of Geary's measure of kurtosis.
}
\usage{
geary(x, na.rm = FALSE)
}
\arguments{
  \item{x}{ a numeric vector, matrix or data frame. }
  \item{na.rm}{ logical. Should missing values be removed? }
}
\details{
  The Geary's kurtosis is computed by dividing average difference between observation
  and the mean by standard deviation of the sample.
}
\references{ 

Geary, R.C. (1936). Moments of the ratio of the mean deviation to the standard deviation
for normal samples. Biometrika, 28, 295-307.

}
\author{ Lukasz Komsta }
\seealso{ \code{\link{kurtosis}}, \code{\link{bonett.test}} }
\examples{
set.seed(1234)
geary(rnorm(1000))
}
\keyword{univar}

