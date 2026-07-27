\name{vix.2013.06.25}
\alias{vix.2013.06.25}
\docType{data}
\title{
VIX Options on 2013-06-25
}
\description{
  This dataset contains VIX options with 57 days to expiration at the end of the business day June 25, 2013.
  On June 25, 2013, VIX closed at 18.21.
}
\usage{data(vix.2013.06.25)}
\format{
  A data frame with 35 observations on the following 13 variables.
  \describe{
    \item{\code{last.c}}{closing call price}
    \item{\code{change.c}}{change in call price from previous day}
    \item{\code{bid.c}}{call bid price}
    \item{\code{ask.c}}{call ask price}
    \item{\code{vol.c}}{call volume}
    \item{\code{openint.c}}{call open interest}
    \item{\code{strike}}{option strike}
    \item{\code{last.p}}{closing put price}
    \item{\code{change.p}}{change in put price from previous day}
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
data(vix.2013.06.25)
}
\keyword{VIX Options}
