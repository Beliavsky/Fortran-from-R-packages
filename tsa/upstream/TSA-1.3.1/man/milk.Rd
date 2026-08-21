\name{milk}
\alias{milk}
\docType{data}
\title{Monthly Milk Production}
\description{ Average monthly milk production per cow in the US,
01/1994 - 12/2005} 
\usage{data(milk)}
\format{
  The format is:
 'ts' int [1:144, 1] 1343 1236 1401 1396 1457 1388 1389 1369 1318 1354 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr "milk"
 - attr(*, "tsp")= num [1:3] 1994 2006   12
}
\examples{
data(milk)
str(milk) 
plot(milk) 
}
\keyword{datasets}
