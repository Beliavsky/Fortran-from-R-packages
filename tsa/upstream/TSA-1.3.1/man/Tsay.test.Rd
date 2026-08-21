\name{Tsay.test}
\alias{Tsay.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Tsay's Test for nonlinearity} 
\description{Carry out Tsay's test for quadratic nonlinearity in a time series. }
\usage{
Tsay.test(x, order, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{time series}
  \item{order}{working linear AR order; if missing, it will be estimated
via the ar function by minimizing AIC}
  \item{\dots}{options to be passed to the ar function}
}
\details{
The null hypothesis is that the true model is an AR process. The AR order, if missing, is estimated by minimizing AIC via the 
ar function, i.e. fitting autoregressive model to the data. The default
fitting method of the ar function is "yule-walker." 
}
\value{
A list containing the following components
  \item{test.stat}{The observed test statistic}
  \item{p.value}{p-value of the test}
  \item{order}{working AR order}
}
\references{Tsay, R. S. (1986), Nonlinearity test for time series,  
Biometrika, 73, 461-466.}
\author{Kung-Sik Chan}
\seealso{\code{\link{Tsay.test}},\code{\link{tlrt}}}
\examples{
data(spots)
Tsay.test(sqrt(spots))
}
\keyword{methods}
