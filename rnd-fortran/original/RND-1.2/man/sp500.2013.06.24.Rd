\name{sp500.2013.06.24}
\alias{sp500.2013.06.24}
\docType{data}
\title{
S&P 500 Index Options on 2013-06-24
}
\description{
  This dataset contains S&P 500 options with 53 days to expiration at the end of the business day June 24, 2013.
  On June 24, 2013, S&P 500 closed at 1573.09.
}
\usage{data(sp500.2013.06.24)}
\format{
  A data frame with 173 observations on the following 9 variables.
  \describe{
    \item{\code{bid.c}}{call bid price}
    \item{\code{ask.c}}{call ask price}
    \item{\code{vol.c}}{call volume}
    \item{\code{openint.c}}{call open interest}
    \item{\code{strike}}{option strike}
    \item{\code{bid.p}}{put bid price}
    \item{\code{ask.p}}{put ask price}
    \item{\code{vol.p}}{put volume}
    \item{\code{openint.p}}{put open interest}
  }
}
\source{
http://www.cboe.com/DelayedQuote/QuoteTableDownload.aspx
}
\examples{
data(sp500.2013.06.24)
}
\keyword{S&P 500}
