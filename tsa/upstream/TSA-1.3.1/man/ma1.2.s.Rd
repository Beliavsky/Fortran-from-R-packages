\name{ma1.2.s}
\alias{ma1.2.s}
\docType{data}
\title{A simulated MA(1) series / time series}
\description{
A simulated MA(1) series with the MA(1) coefficient equal to -0.9.  
}
\usage{data(ma1.2.s)}
\format{
  The format is:
 Time-Series [1:120] from 1 to 120:  1.511  1.821  0.957 -1.538 -2.888 ...
}
\details{
The model is \eqn{Y(t)=e(t)+0.9e(t-1)} where the e's are iid standard normal.  }
\examples{
data(ma1.2.s)
## maybe str(ma1.2.s) ; plot(ma1.2.s) ...
}
\keyword{datasets}
