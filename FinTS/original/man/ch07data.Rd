\name{ch07data}
\docType{data}
\alias{ch07data}
\alias{d.ibm6298wmx}
\alias{d.intc7297}
\title{ Financial time series for Tsay (2005, ch. 7) }
\description{
  Financial time series used in examples in chapter 7.
}
\usage{
data(d.ibm6298wmx)
data(d.intc7297)
}
\format{

  \describe{
    \item{d.ibm6298wmx}{
      a zoo object of 9190 observations on several series relating to
      IBM stock, 1962-07-03 to 1998-12-31:

      \describe{
	\item{dailySimpleRtns}{
	  daily simple returns in percentages of IBM stock
	}
	\item{day}{numbers 1:9190}
	\item{meanCorrectedLogRtns}{mean-corrected log returns}
	\item{Q4}{1 for October, November, December, and 0 otherwise}
	\item{drop2.5pct}{
	  an indicator variable for the behavior of the previous trading
	  day.  Specifically, this is 1 if the meanCorrectedLogRtns for
	  the previous day was at most (-0.025).
	}
	\item{nOfLast5outside2.5pct}{
	  number of the last 5 days for which the meanCorrectedLogRtns
	  exceeded +/-2.5%
	}
	\item{annualTrend}{
	  an annual trend defined as (year-1961)/38.
	}
	\item{GARCH1.1volatility}{
	  a volatility series based on a Gaussian GARCH(1,1) model for
	  the mean-corrected log returns.
	}

	The simpleDailyRtns and the zoo index are from 'd-ibm6298.txt'
	from the book's web site.

	The 'day' and 'meanCorrectedLogRtns' are from 'd-ibmln98wm.txt'.

	The last 5 columns are from 'd-ibml25x.txt';  they are described
	on p. 332 of the book.
      }
    }

    \item{d.intc7297}{
      a zoo object of daily log returns of Intel stock, 1972-12-15 to
      1997-12-31.
    }
  }
}
%\details{
%}
\source{
  \url{https://faculty.chicagobooth.edu/ruey-s-tsay/teaching}
}
\references{
  Ruey Tsay (2005)
  Analysis of Financial Time Series, 2nd ed. (Wiley, ch. 7)
}
\keyword{datasets}
\seealso{
  \code{\link{ch01data}},
  \code{\link{ch02data}},
  \code{\link{ch03data}},
  \code{\link{ch04data}},
  \code{\link{ch05data}},
  \code{\link{ch06data}}
}
%\examples{
%}

