\name{rstandard.Arima}
\alias{rstandard.Arima}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Compute the Standardized Residuals from a Fitted ARIMA Model}
\description{
Computes the internally standardized residuals from a fitted ARIMA model.
}
\usage{
\method{rstandard}{Arima}(model,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{model fitted by the arima function}
  \item{...}{not used; kept here for consistency with the generic method}
}
\details{residuals/(error std. dev.) }
\value{time series of standarized residuals}
\examples{
data(oil.price)
m1.oil=arima(log(oil.price),order=c(0,1,1))
plot(rstandard(m1.oil),ylab='Standardized residuals',type='l')
abline(h=0)
}
\keyword{methods}
