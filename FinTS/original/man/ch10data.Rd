\name{ch10data}
\docType{data}
\alias{ch10data}
\alias{d.hkja}
\alias{m.pfe6503}
\alias{m.mrk6503}
%\alias{m.ibmsp2699}
\alias{d.spcscointc}
\title{ Financial time series for Tsay (2005, ch. 10) }
\description{
  Financial time series used in examples in chapter 10.
}
\usage{
data(d.hkja)
data(m.pfe6503)
data(m.mrk6503)
#data(m.ibmsp2699)
# <- 2 of the 4 columns in m.ibmsp2699ln
# documented with ch08data
data(d.spcscointc)
}
\format{
  One data set used in chapter 10 is also used earlier:
  'm.ibmsp2699' is the first 2 of the 4 columns of 'm.ibmsp2699ln' used
  in chapter 8.

  The other data sets used in chapter 10 are as follows:

  \describe{
    \item{d.hkja}{
      zoo object giving the daily log returns of HK and Japan market
      indices from 1996-01-01 through 1997-05-05 (used in Example 10.1).
    }
    \item{m.pfe6503, m.mrk6503}{
      zoo objects giving the monthly simple returns including dividends
      of Pfizer and Merk stocks.
    }
    \item{d.spcscointc}{
      data.frame giving 2275 daily log returns of three items from
      January 2, 1991 through December 31, 1999:

      \describe{
	\item{SP500}{Standard & Poor's 500 index}
	\item{Cisco}{Cisco stock}
	\item{Intel}{Intel stock}
      }

      NOTE:  This date range seems to include 2280 trading days in the
      New York Stock Exchange.  Since the file on the book's web site
      did not include dates and since there appear to be more trading
      days than observations, dates are not currently provided with
      these observations.  This may change with a future revision of
      this package.
    }
  }
}
%\details{}
\source{
  \url{https://faculty.chicagobooth.edu/ruey-s-tsay/teaching}
}
\references{
  Ruey Tsay (2005) Analysis of Financial Time Series, 2nd ed. (Wiley,
  ch. 10)
}
\keyword{datasets}
\seealso{
  \code{\link{ch01data}},
  \code{\link{ch02data}},
  \code{\link{ch03data}},
  \code{\link{ch04data}},
  \code{\link{ch05data}},
  \code{\link{ch06data}},
  \code{\link{ch07data}},
  \code{\link{ch08data}},
  \code{\link{ch09data}}
}
%\examples{
%}
