\name{ch09data}
\docType{data}
\alias{ch09data}
\alias{m.fac9003}
\alias{m.cpice16.dp7503}
\alias{m.barra.9003}
\alias{m.5cln}
%\alias{m.bnd}
\alias{m.apca0103}
\title{ Financial time series for Tsay (2005, ch. 9) }
\description{
  Financial time series used in examples in chapter 9.
}
\usage{
data(m.fac9003)
data(m.cpice16.dp7503)
data(m.barra.9003)
data(m.5cln)
#data(m.bnd) <- documented with ch08, also used in ch09
data(m.apca0103)
}
\format{

  \describe{
    \item{m.fac9003}{
      a zoo object of 168 observations giving simple excess returns of
      13 stocks and the Standard and Poor's 500 index over the monthly
      series of three-month Treasury bill rates of the secondary market
      as the risk-free rate from January 1990 to December 2003.  (These
      numbers are used in Table 9.1.)

      \describe{
	\item{AA}{Alcoa}
	\item{AGE}{A. G. Edwards}
	\item{CAT}{Caterpillar}
	\item{F}{Ford Motor}
	\item{FDX}{FedEx}
	\item{GM}{General Motors}
	\item{HPQ}{Hewlett-Packard}
	\item{KMB}{Kimberly-Clark}
	\item{MEL}{Mellon Financial}
	\item{NYT}{New York Times}
	\item{PG}{Proctor & Gamble}
	\item{TRB}{Chicago Tribune}
	\item{TXN}{Texas Instruments}
	\item{SP5}{Standard & Poor's 500 index}
      }
    }

    \item{m.cpice16.dp7503}{
      a zoo object of 168 monthly on two macroeconomic variables from
      January 1975 through December 2002 (p. 412):

      \describe{
	\item{CPI}{
	  consumer price index for all urban consumers:  all items and
	  with index 1982-1984 = 100
	}
	\item{CE16}{
	  Civilian employment numbers 16 years and over:  measured in
	  thousands
	}
      }
    }
    \item{m.barra.9003}{
      a zoo object giving monthly excess returns of ten stocks from
      January 1990 through December 2003:

      \describe{
	\item{AGE}{A. G. Edwards}
	\item{C}{Citigroup}
	\item{MWD}{Morgan Stanley}
	\item{MER}{Merrill Lynch}
	\item{DELL}{Dell, Inc.}
	\item{IBM}{International Business Machines}
	\item{AA}{Alcoa}
	\item{CAT}{Caterpillar}
	\item{PG}{Proctor & Gamble}
      }
    }
    \item{m.5cln}{
      a zoo object giving monthly log returns in percentages of 5 stocks
      from January 1990 through December 1999:

      \describe{
	\item{IBM}{International Business Machines}
	\item{HPQ}{Hewlett-Packard}
	\item{INTC}{Intel}
	\item{MER}{Merrill Lynch}
	\item{MWD}{Morgan Stanley Dean Witter}
      }
    }
    \item{m.apca0103}{
      data.frame of monthly simple returns of 40 stocks from January
      2001 through December 2003, discussed in sect. 9.6.2, pp. 437ff.

      \describe{
	\item{CompanyID}{5-digit company identification code}
	\item{date}{the last workday of the month}
	\item{return}{in percent}
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
\examples{
data(m.apca0103)
dim(m.apca0103)
# 1440 3;  1440 = 40*36
# Are the dates all the same?
sameDates <- rep(NA, 39)
for(i in 1:39)
    sameDates[i] <- with(m.apca0103,
                         all.equal(date[1:36], date[(i*36)+1:36]))
stopifnot(all(sameDates))
M.apca0103 <- with(m.apca0103, array(return, dim = c(36, 40), dimnames =
    list(NULL, paste("Co", CompanyID[seq(1, 1440, 36)], sep=""))))
}
