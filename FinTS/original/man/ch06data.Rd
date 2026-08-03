\name{ch06data}
\docType{data}
\alias{ch06data}
\alias{d.ibmy98}
\alias{d.cscoy99}
\title{ Financial time series for Tsay (2005, ch. 6) }
\description{
  Financial time series used in examples in chapter 6.
}
\usage{
data(d.ibmy98)
data(d.cscoy99)
}
\format{
  Objects of class zoo giving returns for each trading day for different
  periods

  \describe{
    \item{d.ibmy98}{
      Zoo object giving daily  simple returns of IBM stock for each
      trading day in 1998.
    }
    \item{d.cscoy99}{
      Zoo object giving daily log returns of Cisco stock for each
      trading day in 1999.
    }
  }
}
%\details{}
\source{
  \url{https://faculty.chicagobooth.edu/ruey-s-tsay/teaching}
}
\references{
  Ruey Tsay (2005)
  Analysis of Financial Time Series, 2nd ed. (Wiley, ch. 6)
}
\keyword{datasets}
\seealso{
  \code{\link{ch01data}},
  \code{\link{ch02data}},
  \code{\link{ch03data}},
  \code{\link{ch04data}},
  \code{\link{ch05data}}
}
%\examples{}
