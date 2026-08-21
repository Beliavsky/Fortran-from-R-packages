\name{LB.test}
\alias{LB.test}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Portmanteau Tests for Fitted ARIMA models}
\description{
This function modifies the Box.test function in the stats package, and it
computes the Ljung-Box or Box-Pierce tests checking whether or not the
residuals appear to be white noise. 
}
\usage{
LB.test(model, lag = 12, type = c("Ljung-Box", "Box-Pierce"), no.error = FALSE,
 omit.initial = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{model fit from the arima function}
  \item{lag}{number of lags of the autocorrelation of the residuals to be
included in the test statistic. (default=12)}
  \item{type}{either Ljung-Box or Box-Pierce}
  \item{no.error}{a system variable; normally it is not changed}
  \item{omit.initial}{if true, (d+Ds) initial residuals are omitted from 
the test}
}
\value{
a list:
  \item{statistics}{test statistic}
  \item{p.value}{p-value}
\item{parameter}{d.f. of the Chi-square test}
\item{lag}{no of lags}
}
\author{Kung-Sik Chan, based on A. Trapletti's work  on the Box.test function in the stats package}
\references{
Box, G. E. P. and Pierce, D. A. (1970), Distribution of residual correlations in autoregressive-integrated moving average time series models. Journal of the American Statistical Association, 65, 15091526. 

Ljung, G. M. and Box, G. E. P. (1978), On a measure of lack of fit in time series models. Biometrika 65, 553564. 

}
\examples{
data(color)
m1.color=arima(color,order=c(1,0,0))
LB.test(m1.color)
}
\keyword{methods}
