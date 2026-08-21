\name{fitted.Arima}
\alias{fitted.Arima}
\title{Fitted values of an arima model.}
\description{
Computes the fitted values of an arima model.}
\usage{
\method{fitted}{Arima}(object,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{a fitted model from the arima function.}
  \item{...}{other arguments; not used here but kept to be consistent with the 
generic method}
}
\value{
fitted values}
\author{Kung-Sik Chan}
\seealso{ \code{\link{arima}}}
\examples{
data(hare)
hare.m1=arima(sqrt(hare),order=c(3,0,0))
fitted(hare.m1)
}
\keyword{methods}
