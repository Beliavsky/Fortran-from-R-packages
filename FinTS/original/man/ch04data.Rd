\name{ch04data}
\docType{data}
\alias{ch04data}
\alias{m.unrate}
%\alias{d.ibmvwewsp6203}
%\alias{m.3m4697}
%\alias{q.gnp4791}
\alias{w.3mtbs7097}
\alias{m.ibmln2699}
\alias{q.unemrate}
\title{ Financial time series for Tsay (2005, ch. 4) }
\description{
  Financial time series used in examples in chapter 4.
}
\usage{
data(m.unrate)
#d.ibmvwewsp6203
#m.3m4697
#q.gnp4791
data(w.3mtbs7097)
data(m.ibmln2699)
data(q.unemrate)
}
\format{
  Three data stes used in chapter 4 are also used earlier:
  'd.ibmvwewsp6203' is used in chapter 1, and 'm.3m4697' and 'q.gnp4791'
  are used in chapter 2;  these three data objects are docuemented in
  'ch01data' or 'ch02data'.

  The other data sets used in chapter 4 are as follows:

  \describe{
    \item{m.unrate}{
      zoo object giving the monthly US civilian unemployment rate from
      1948 through 2004.
    }
    \item{w.3mtbs7097}{
      zoo object giving the US weekly 3-month treasury bill rate in the
      secondary market from 1970 through 1997.
    }
    \item{m.ibmln2699}{
      zoo object giving the monthly log returns in percentages of IBM
      stock from 1926 through 1999.
    }
    \item{q.unemrate}{
      zoo object giving the US quarterly unemployment rate seasonally
      adjusted from 1948 through 1993.
    }
  }
}
%\details{
%}
\source{
  \url{https://faculty.chicagobooth.edu/ruey-s-tsay/teaching}
}
\references{
  Ruey Tsay (2005) Analysis of Financial Time Series, 2nd ed. (Wiley,
  ch. 4)
}
\keyword{datasets}
\seealso{
  \code{\link{ch01data}}
  \code{\link{ch02data}}
}
%\examples{
%}

