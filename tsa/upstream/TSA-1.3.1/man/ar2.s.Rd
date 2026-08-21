\name{ar2.s}
\alias{ar2.s}
\docType{data}
\title{Asimulated AR(2) series / time series}
\description{
Asimulated AR(2) series with AR coefficients being equal to 1.5 and -0.75} 
\usage{data(ar2.s)}
\format{
  The format is:
 Time-Series [1:120] from 1 to 120: -2.064 -1.937  0.406  2.039  2.953 ...
}
\details{
The model is Y(t)=1.5*Y(t-1)-0.75*Y(t-2)+e(t) where the e's are iid
standard normal random variables.
}
\examples{
data(ar2.s)
## maybe str(ar2.s) ; plot(ar2.s) ...
}
\keyword{datasets}
