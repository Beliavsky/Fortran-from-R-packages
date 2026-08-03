\name{ch11data}
\docType{data}
\alias{ch11data}
\alias{aa.3rv}
\title{ Financial time series for Tsay (2005, ch. 11) }
\description{
  Financial time series used in examples in chapter 11.
}
\usage{
data(aa.3rv)
# m.fac9003 described in ch09data
# q.jnj described in dh02data
}
\format{
  The text of chapter 11 considers one data set not used in previous
  chapters plus two that are.  Monthly excess returns of GM stock are
  used in Table 9.1 of Chapter 9.  Quarterly earnings of Johnson and
  Johnson are used in Chapter 2.

  The data set introduced with Chapter 11 is as follows:

  \describe{
    \item{aa.3rv.}{
      a zoo object of daily 5, 10, and 20 minute realized volatility of
      Alcoa stock from 2003-01-02 through 2004-05-07.
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
  ch. 11)
}
\keyword{datasets}
\seealso{
  \code{\link{ch09data}},
  \code{\link{ch02data}}
}
%\examples{
%}
