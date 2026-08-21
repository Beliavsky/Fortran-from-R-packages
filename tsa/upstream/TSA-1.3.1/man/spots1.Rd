\name{spots1}
\alias{spots1}
\docType{data}
\title{
Annual international sunspot numbers
}
\description{
Annual international sunspot numbers, NOAA National Geophysical
Data Center, 1700 - 2005.
}
\usage{data(spots1)}
\format{
  The format is:
 ts [1:306, 1] 5 11 16 23 36 58 29 20 10 8 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr "spots"
 - attr(*, "tsp")= num [1:3] 1700 2005    1
}
\source{
\code{ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SUNSPOT_NUMBERS/YEARLY.PLT}
}
\examples{
data(spots1)
## maybe str(spots1) ; plot(spots1) ...
}
\keyword{datasets}
