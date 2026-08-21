\name{ma2.s}
\alias{ma2.s}
\docType{data}
\title{A simulated MA(2) series}
\description{
A simulated MA(2) series with MA coefficients being 1 and -0.6.
}
\usage{data(ma2.s)}
\format{
  The format is:
 Time-Series [1:120] from 1 to 120: -0.4675  0.0815  0.9938 -2.6959  2.8116 ...
}
\details{
The model is \eqn{Y(t)=e(t)-e(t-1)+0.6*e(t-2)} where the e's are iid standard normal
random variables. 
}
\examples{
data(ma2.s)
## maybe str(ma2.s) ; plot(ma2) ...
}
\keyword{datasets}
