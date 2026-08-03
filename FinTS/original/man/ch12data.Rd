\name{ch12data}
\docType{data}
\alias{ch12data}
\alias{w.gs3n1c}
\alias{w.gs3c}
\alias{m.sp6299}
\alias{m.ibmspln6299}
\alias{m.sp5.6204}
\alias{m.geln}
\title{ Financial time series for Tsay (2005, ch. 12) }
\description{
  Financial time series used in examples in chapter 12.
}
\usage{
data(w.gs3n1c)
data(w.gs3c)
data(m.sp6299)
data(m.ibmspln6299)
data(m.sp5.6204)
data(m.geln)
}
\format{

  \describe{
    \item{w.gs3n1c}{
      a zoo object of the change series of weekly US interest rates (3
      and 1 year maturities) from  Jan. 5, 1962, to Sep. 10, 1999.  This
      was obtained via
      
      \code{diff(window(w.gs1n36299, start=as.Date("1962-01-05"), end=as.Date("1999-09-10")))[, 2:1]}

      to get the dates with the data. Then 'all.equal' confirmed that
      these numbers matched those in the file read from the web site
      (which did not have dates).

      These are used in Example 12.1, pp. 556ff.
    }
    \item{w.gs3c}{
      
      a zoo object giving the change series of weekly US 3-year maturity
      interest rates from March 18, 1988, to Sept. 10, 1999.  This was
      obtained via

      \code{window(w.gs3n1c[, 1], start=as.Date("1988-03-18"), end = as.Date("1999-09-10"))}.

      Then 'all.equal' confirmed that these numbers matched those read
      from the web site.

      These data are used in Example 12.2, pp. 564ff.
    }
    \item{m.sp6299}{
      Monthly log returns of S&P 500 index from January 1962 to December
      1999.  These data are used in Example 12.3, pp. 569ff.

      These data are a subset of 'm.ibmspln', used in chapter 3.  That
      series has dates, which were not provided in the file associated
      with this series on the book's web site.  Moreover, the file with
      chapter 12 has only 4 significant digits where the earlier file
      has 6.  Since the other data are otherwise identical, this
      'm.sp6299' was constructed as

      \code{window(m.ibmspln[, 2], start = yearmon(1962), end = yearmon(1999+11/12))}.
    }
    \item{m.ibmspln6299}{
      Monthly log returns of IBM stock and the S&P 500 index from
      January 1962 to December 1999.  These data are used in Example
      12.4, pp. 573ff.

      These data are an expansion of 'm.sp6299' and were similarly
      obtained from 'm.ibmspln'.
    }
    \item{m.sp5.6204}{
      Monthly log returns of S&P 500 index from January 1962 to November
      1999.  These data are used in Example 12.5, pp. 586ff.
    }
    \item{m.geln}{
      Monthly log returns of GE stock from January 1926 to December
      1999.  These data are used in Example 12.6, pp. 591ff.
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
  ch. 12)
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
  \code{\link{ch09data}},
  \code{\link{ch10data}},
  \code{\link{ch11data}}
}
%\examples{
%}
