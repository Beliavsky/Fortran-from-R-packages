\name{ar1.2.s}
\alias{ar1.2.s}
\docType{data}
\title{A simulated AR(1) series}
\description{
A simulated AR(1) series with the AR coefficient equal  to 0.4.

}
\usage{data(ar1.2.s)}
\format{
  The format is:
Time-Series [1:60] from 1 to 60: -0.0678  1.4994  0.4888  0.3987 -0.5162 ...
}
\details{
The model is Y(t)=0.4*Y(t-1)+e(t) where the e's are iid standard normal.
}
\examples{
data(ar1.2.s)
## maybe str(ar1.2.s) ; plot(ar1.2.s) ...
}
\keyword{datasets}
