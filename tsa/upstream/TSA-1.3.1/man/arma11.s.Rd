\name{arma11.s}
\alias{arma11.s}
\docType{data}
\title{ A Simulated ARMA(1,1) Series/ time series}
\description{
A simulated ARMA(1,1) series with the model given by:
\eqn{y_t=0.6*y_{t-1}+e_t+0.3*e_{t-1}} where the e's are iid standard normal 
random variables.
}
\usage{data(arma11.s)}
\format{
  The format is:
 Time-Series [1:100] from 1 to 100: -0.765  1.297  0.668 -1.607 -0.626 ...
}
\examples{
data(arma11.s)
## maybe str(arma11.s) ; plot(arma11.s) ...
}
\keyword{datasets}
