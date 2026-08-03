\name{ch08data}
\docType{data}
\alias{ch08data}
\alias{m.ibmsp2699ln}
\alias{m.bnd}
\alias{m.gs1n3.5301}
\alias{w.tb3n6ms}
\alias{sp5may}
\title{ Financial time series for Tsay (2005, ch. 8) }
\description{
  Financial time series used in examples in chapter 8.
}
\usage{
data(m.ibmsp2699ln)
data(m.bnd)
data(m.gs1n3.5301)
data(w.tb3n6ms)
data(sp5may)
}
\format{

  \describe{
    \item{m.ibmsp2699ln}{
      zoo object giving monthly simple and log returns of IBM stock and
      the Standard and Poor's 500 from 1926 through 1999.  (This
      combines files 'm-ibmsp2699.txt' and 'm-ibmspln.txt' from the
      book's web site.)
    }
    \item{m.bnd}{
      zoo object giving the monthly simple returns of 30, 20, 10, 5 and
      1 year maturity bonds from 1942 through 1999.
    }
    \item{m.gs1n3.5301}{
      zoo object giving 1 and 3 year US Treasury constant maturity
      interest rates from April 1953 to January 2001 (used in Example
      8.6, pp. 373ff).
    }
    \item{w.tb3n6ms}{
      zoo object giving weekly 3 and 6 month US Treasury Bill interest
      rates from 1958-12-12 to 2004-08-06 (used in Sect. 8.6.5,
      pp. 385ff).
    }
    \item{sp5may}{
      A data.frame of 7061 observations on 4 variables based on
      minute-by-minute observations of the Standard and Poor's 500
      Futures and prices in May 1993.

      These data are used, after some processing, in Tsay(Sect. 8.7.2,
      pp. 392ff).  Unfortunately, it's not yet clear what these numbers
      are.  The following is a current guess and will doubtless change
      in the future.

      \describe{
	\item{logFuture}{
	  logarithms of June Futures contracts traded at the Chicago
	  Mercantile Exchange.  The first difference of this series
	  appears to be plotted in Figure 8.16(a), after replacing '10
	  extreme values (5 on each side) by the simple average of their
	  two nearest neighbors.' (p. 392)
	}
	\item{logPrice}{
	  logarithms of Standard and Poor's 500 price levels.  The first
	  differences of this series appears to be plotted in Figure
	  8.16(b), after adjustment similar to that for 'logFuture'.
	}
	\item{dailyAvgSomething}{
	  numbers that assume 19 distinct levels separated by 18
	  discrete jumps.  The name of this will likely change whenever
	  more information about it can be obtained for this
	  documentation.
	}
	\item{day}{
	  index for the 19 distinct levels assumed by
	  'dailyAvgSomething'.  This is probably the trading day in May
	  1993.  However, there appear to have been 20 trading days in
	  that month, so if these 19 levels do correspond to trading
	  days, it's not clear which date is missing.
	}

	These data were analyzed by Forbes, Kalb, and Kofman (1999);
	Tsay (1998) was also referenced with the discussion of the
	analysis of these data.
      }
    }
  }

}
%\details{
%}
\source{
  \url{https://faculty.chicagobooth.edu/ruey-s-tsay/teaching}
}
\references{
  Forbes, C. S., Kalb, G. R. J., and Kofman, P. (1999) 'Bayesian
  Arbitrage Threshold Analysis', Journal of Business and Economic
  Statistics, 17:  364-372.

  Ruey Tsay (1998) 'Testing and Modeling Multivariate Threshold Models',
  Journal of the American Statistical Association, 93:  1188-1202.

  Ruey Tsay (2005) Analysis of Financial Time Series, 2nd ed. (Wiley,
  ch. 8)
}
\keyword{datasets}
\seealso{
  \code{\link{ch01data}},
  \code{\link{ch02data}}
}
%\examples{
%}
