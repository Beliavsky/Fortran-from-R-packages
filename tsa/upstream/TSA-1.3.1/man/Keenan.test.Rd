\name{Keenan.test}
\alias{Keenan.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Keenan's one-degree test for nonlinearity}
\description{
Carry out Keenan's 1-degree test for nonlinearity against the null
hypothesis that the time series follows some AR process.}
\usage{
Keenan.test(x, order, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{time series}
  \item{order}{working AR order; if missing, it is estimated by minimizing
AIC via the ar function.}
  \item{\dots}{user-supplied options to the ar function.}
}
\details{
The test is designed to have  optimal local power against depature 
from the linear autoregressive function in the direction of the square of
the linear autoregressive function.
}
\value{
A list containing the following components
  \item{test.stat}{The observed test statistic}
  \item{p.value}{p-value of the test}
  \item{order}{working AR order}
}
\references{Keenan, D. M. (1985), A Tukey nonadditivity-type test for time series Nonlinearity,  Biometrika,  72,   39-44.}
\author{Kung-Sik Chan}
\seealso{\code{\link{Tsay.test}},\code{\link{tlrt}}}
\examples{
data(spots)
Keenan.test(sqrt(spots))
}
\keyword{methods}
