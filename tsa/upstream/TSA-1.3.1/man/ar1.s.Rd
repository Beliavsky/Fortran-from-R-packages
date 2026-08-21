\name{ar1.s}
\alias{ar1.s}
\docType{data}
\title{A simulated AR(1) series}
\description{
A simulated AR(1) series with the AR coefficient equal  to 0.9.

}
\usage{data(ar1.s)}
\format{
  The format is:
 Time-Series [1:60] from 1 to 60: -1.889 -1.691 -1.962 -0.566 -0.627 ...
}
\details{
The model is Y(t)=0.9*Y(t-1)+e(t) where the e's are iid standard normal.
}
\examples{
data(ar1.s)
## maybe str(ar1.s) ; plot(ar1.s) ...
}
\keyword{datasets}
