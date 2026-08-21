\name{airmiles}
\alias{airmiles}
\docType{data}
\title{Monthly Airline Passenger-Miles in the US}
\description{
Monthly U.S. airline passenger-miles: 01/1996 - 05/2005.}
\usage{data(airmiles)}
\format{
  The format is:
 'ts' int [1:113, 1] 30983174 32147663 38342975 35969113 36474391 38772238 40395657 41738499 33580773 36389842 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr "airmiles"
 - attr(*, "tsp")= num [1:3] 1996 2005   12
}
\source{\code{www.bts.gov/xml/air_traffic/src/index.xml#MonthlySystem}}
\examples{
data(airmiles)
## maybe str(airmiles) ; plot(airmiles) ...
}
\keyword{datasets}
