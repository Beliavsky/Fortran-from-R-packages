\name{bluebird}
\alias{bluebird}
\docType{data}
\title{Blue Bird Potato Chip Data}
\description{
Weekly unit sales (log-transformed) of Bluebird standard potato chips 
(New Zealand)
and their price for 104 weeks.}
\usage{data(bluebird)}
\format{
  The format is:
 mts [1:104, 1:2] 11.5 11.5 11.8 11.9 11.3 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:2] "log.sales" "price"
 - attr(*, "tsp")= num [1:3] 1 104 1
 - attr(*, "class")= chr [1:2] "mts" "ts"
}
\source{www.stat.auckland.ac.nz/~balemi/Assn3.xls}
\examples{
data(bluebird)
## maybe str(bluebird) ; plot(bluebird) ...
}
\keyword{datasets}
