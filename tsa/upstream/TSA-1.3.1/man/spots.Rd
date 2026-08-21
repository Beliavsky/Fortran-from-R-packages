\name{spots}
\alias{spots}
\docType{data}
\title{Relative annual sunspot number / time series}
\description{ Annual American (relative) sunspot numbers collected
from 1945 to 2007.
 The annual (relative) sunspot number is a weighted average of
solar activities measured from a network of observatories.
}
\usage{data(spots)}
\format{
  The format is:
 Time-Series [1:61] from 1945 to 2005:  32.3  99.9 170.9 166.6 174.1 ...
}
\source{
\code{http://www.ngdc.noaa.gov/stp/SOLAR/ftpsunspotnumber.html#american}
}
\references{ "Time Series Analysis, with Applications in R" by J.D. Cryer and K.S. Chan}
\examples{
data(spots)
## maybe str(spots) ; plot(spots) ...
}
\keyword{datasets}
