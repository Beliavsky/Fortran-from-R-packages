\name{sp500.2013.04.19}
\alias{sp500.2013.04.19}
\docType{data}
\title{
S&P 500 Index Options on 2013-04-19
}
\description{
  This dataset contains S&P 500 options with 62 days to expiration at the end of the business day April 19, 2013.  
  On April 19, 2013, S&P 500 closed at 1555.25.
}
\usage{data(sp500.2013.04.19)}
\format{
  A data frame with 171 observations on the following 19 variables.
  \describe{
    \item{\code{bidsize.c}}{call bid size}
    \item{\code{bid.c}}{call bid price}
    \item{\code{ask.c}}{call ask price}
    \item{\code{asksize.c}}{call ask size}
    \item{\code{chg.c}}{change in call price}
    \item{\code{impvol.c}}{call implied volatility}
    \item{\code{vol.c}}{call volume}
    \item{\code{openint.c}}{call open interest}
    \item{\code{delta.c}}{call delta}
    \item{\code{strike}}{option strike}
    \item{\code{bidsize.p}}{put bid size}
    \item{\code{bid.p}}{put bid price}
    \item{\code{ask.p}}{put ask price}
    \item{\code{asksize.p}}{put ask size}
    \item{\code{chg.p}}{change in put price}
    \item{\code{impvol.p}}{put implied volatility}
    \item{\code{vol.p}}{put volume}
    \item{\code{openint.p}}{put open interest}
    \item{\code{delta.p}}{put delta}
  }
}
\source{
http://www.cboe.com/DelayedQuote/QuoteTableDownload.aspx
}
\examples{
data(sp500.2013.04.19)
}
\keyword{S&P 500}
