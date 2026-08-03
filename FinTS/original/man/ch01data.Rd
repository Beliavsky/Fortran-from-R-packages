\name{ch01data}
\docType{data}
\alias{ch01data}
\alias{Table1.2}
\alias{d.ibmvwewsp6203}
\alias{d.intc7303}
\alias{d.3m6203}
\alias{d.msft8603}
\alias{d.c8603}
\alias{m.ibmvwewsp2603}
\alias{m.intc7303}
\alias{m.3m4603}
\alias{m.msft8603}
\alias{m.c8603}
\alias{m.gs10}
\alias{m.gs1}
\alias{d.fxjp00}
\alias{m.fama.bond5203}
\alias{m.gs3}
\alias{m.gs5}
\alias{w.tb3ms}
\alias{w.tb6ms}
\title{ Financial time series for Tsay (2005, ch. 1) }
\description{
  Financial time series used in examples in chapter 1.
}
\usage{
data(d.ibmvwewsp6203)
data(d.intc7303)
data(d.3m6203)
data(d.msft8603)
data(d.c8603)
data(m.ibmvwewsp2603)
data(m.intc7303)
data(m.3m4603)
data(m.msft8603)
data(m.c8603)
data(m.gs10)
data(m.gs1)
data(d.fxjp00)
data(m.fama.bond5203)
data(m.gs3)
data(m.gs5)
data(w.tb3ms)
data(w.tb6ms)
}
\format{
  Objects of class zoo giving simple returns for each trading period
  (day, week or month) for different periods, with different start dates
  but typically running to the end of 2003.

  \describe{
    \item{d.ibmvwewsp6203, m.ibmvwewsp2603}{
      Zoo objects with 4 columns (IBM, VW, EW, and SP).  Daily data
      starts with 1962-07-03.  Monthly data starts with 1926-01-30.
    }
    \item{d.intc7303, m.intc7303}{
      Matrices of class zoo with a single column "Intel" starting from
      January 1973.
    }
    \item{d.3m6203, m.3m6203}{
      Matrices of class zoo with a single column "MMM".  Daily data
      starts with 1962-07-03.  Monthly data starts with 1946-02-28.
    }
    \item{d.msft8603, m.msft8603}{
      Matrices of class zoo with a single column "MSFT" starting from
      1906-03-14.
    }
    \item{d.c8603, m.c8603}{
      Matrix of class zoo with a single column "C" starting from
      1986-10-30.
    }
    \item{m.gs10, m.gs1}{
      Monthly 10-yr  and 1-yr Treasury constant maturity rates
      (4/53-3/04)
    }
    \item{d.fxjp00}{
      Daily exchange rate between U.S. dollar and Japanese yen
    }
    \item{m.fama.bond5203}{
      Monthly bond returns as follows:
      \describe{
	\item{m1.12}{1-12m}
	\item{m24.36}{24-36m}
	\item{m48.60}{48-60m}
	\item{m61.120}{61-120m}
      }
    }
    \item{m.gs3, m.gs5}{
      Monthly 3-yr and 5-yr Treasury constant maturity rates
    }
    \item{w.tb3ms, w.tb6ms}{
      Weekly Treasury Bill rates
    }
  }
}
\details{
  The first 16 of these objects contain daily and monthly simple returns
  for 8 financial time series analyzed Tsay (2005, Table1.2).  These 8
  are SP (Standard & Poors), EW, IBM, Intel, Microsoft, and Citi-Group,
  beginning at different times and running to the end of 2003.

  The others are used elsewhere in chapter 1.
}
\source{
  \url{https://faculty.chicagobooth.edu/ruey-s-tsay/teaching}
}
\references{
  Ruey Tsay (2005) Analysis of Financial Time Series, 2nd ed. (Wiley,
  ch. 1)
}
\keyword{datasets}
\seealso{ \code{\link{FinTS.stats}} }
\examples{
# First half of Table 1.2:
data(d.ibmvwewsp6203)
data(d.intc7303)
data(d.3m6203)
data(d.msft8603)
data(d.c8603)
(Daily.Simple.Returns.pct <- rbind(
    SP = FinTS.stats(100*d.ibmvwewsp6203[, "SP"]),
    VW = FinTS.stats(100*d.ibmvwewsp6203[, "VW"]),
    EW = FinTS.stats(100*d.ibmvwewsp6203[, "EW"]),
    IBM= FinTS.stats(100*d.ibmvwewsp6203[, "IBM"]),
    Intel=FinTS.stats(100*d.intc7303[, "Intel"]),
    MMM= FinTS.stats(100*d.3m6203[, "MMM"]),
    MSFT=FinTS.stats(100*d.msft8603[, 'MSFT']),
    C  = FinTS.stats(100*d.c8603[, "C"])
) )

(Daily.log.Returns.pct <- rbind(
    SP = FinTS.stats(100*log(1+d.ibmvwewsp6203[, "SP"])),
    VW = FinTS.stats(100*log(1+d.ibmvwewsp6203[, "VW"])),
    EW = FinTS.stats(100*log(1+d.ibmvwewsp6203[, "EW"])),
    IBM= FinTS.stats(100*log(1+d.ibmvwewsp6203[, "IBM"])),
    Intel=FinTS.stats(100*log(1+d.intc7303[,"Intel"])),
    MMM= FinTS.stats(100*log(1+d.3m6203[, "MMM"])),
    MSFT=FinTS.stats(100*log(1+d.msft8603[, 'MSFT'])),
    C  = FinTS.stats(100*log(1+d.c8603[, "C"]))
) )
}
