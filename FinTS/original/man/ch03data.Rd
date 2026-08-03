\name{ch03data}
\docType{data}
\alias{ch03data}
%\alias{m.intc7303}
\alias{exch.perc}
\alias{sp500}
%\alias{m.ibm2697}
%\alias{d.ibmvwewsp6203}
\alias{m.ibmspln}
\alias{m.ibmsplnsu}
\alias{d.sp8099}
\title{ Financial time series for Tsay (2005, ch. 3) }
\description{
  Financial time series used in examples in chapter 3.
}
\usage{
#m.intc7303
data(exch.perc)
data(sp500)
#m.ibm2697
#d.ibmvwewsp6203
data(m.ibmspln)
data(m.ibmsplnsu)
data(d.sp8099)
}
\format{
  Three data stes used in chapter 3 are also used in chapter 1 or 2 and
  are documented with 'ch01data' or 'ch02data':  In particular,
  'm.intc7303' and 'd.ibmvwewsp6203' are used in chapters 1 and 3 and
  are documented with 'ch01data';  'm.ibm2697' is used in chapters 2 and
  3 is documented with ch02data.

  The other data sets used in chapter 3 are as follows:

  \describe{
    \item{exch.perc}{
      numeric vector of length 2497 giving percentage changes in the
      exchange rate between the German mark and the US dollar in 10
      minute intervals, June 5-19, 1989.  (The book describes analyses
      of 2488 observations.  If these 2497 observations are plotted, it
      is difficult to see any differences from Figure 3.2.)
    }
    \item{sp500}{
      object of class 'zooreg' giving the monthly excess returns of the
      S&P 500 index starting from 1926.  This zooreg object is labeled
      assuming it starts in January, though the book does not say
      whether it starts in January or just some time in 1926.  (Many of
      the files included date with the data, but 'sp500.dat' did not.)
    }
    \item{m.ibmspln}{
      object of class 'zooreg' giving the monthly log returns of IBM
      stock and S&P 500 index from January 1926 to December 1999 for 888
      observations.  NOTE:  The examples in the book use only the first
      864 of these observations.
    }
    \item{m.ibmsplnsu}{
      same as 'm.ibmspln' but with a third column 'summer' that is 1 in
      June, July and August, and 0 otherwise.
    }
    \item{d.sp8099}{
      zoo object giving the average daily returns of the S&P 500 from
      1980 through 1999.
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
  ch. 3)
}
\keyword{datasets}
\seealso{
  \code{\link{ch01data}},
  \code{\link{ch02data}}
}
%\examples{
%}

