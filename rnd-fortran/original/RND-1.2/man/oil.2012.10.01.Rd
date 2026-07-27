\name{oil.2012.10.01}
\alias{oil.2012.10.01}
\docType{data}
\title{
West Texas Intermediate Crude Oil Options on 2013-10-01
}
\description{
  This dataset contains West Texas Intermediate (WTI) crude oil options with 43 days to expiration at the end of 
  the business day October 1, 2012.  On October 1, 2012, WTI closed at 92.44.
}
\usage{data(oil.2012.10.01)}
\format{
  A data frame with 332 observations on the following 7 variables.
  \describe{
    \item{\code{type}}{a factor with levels \code{C} for call option \code{P} for put option}
    \item{\code{strike}}{option strike}
    \item{\code{settlement}}{option settlement price}
    \item{\code{openint}}{option open interest}
    \item{\code{volume}}{trading volume}
    \item{\code{delta}}{option delta}
    \item{\code{impliedvolatility}}{option implied volatility}
  }
}

\source{
CME posts sample data at:
http://www.cmegroup.com/market-data/datamine-historical-data/endofday.html
}
\examples{
data(oil.2012.10.01)
}
\keyword{WTI Crude oil options}
