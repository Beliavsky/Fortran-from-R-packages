\name{veilleux}
\alias{veilleux}
\docType{data}
\title{An experimental prey-predator time series}
\description{
A data frame consisting of bivariate time series from an experiment for 
studying prey-predator dynamics. The first time series consists of
  the numbers of
prey individuals (Didinium natsutum) per ml measured
every twelve hours over a period of 35 days; the second time series consists
of the corresponding number of predators (Paramecium aurelia) per ml.
}
\usage{data(veilleux)}
\format{
  The format is:
 mts [1:71, 1:2]  15.7  53.6  73.3  93.9 115.4 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr [1:2] "Didinium" "Paramecium"
 - attr(*, "tsp")= num [1:3] 0 35 2
 - attr(*, "class")= chr [1:2] "mts" "ts"
}
\source{\code{http://www-personal.buseco.monash.edu.au/~hyndman/TSDL/}
}
\references{
Veilleux (1976) "The analysis of a predatory interaction between 
Didinium and Paramecium", Masters thesis, University of Alberta. 

Jost & Ellner (2000) "Testing for predator dependence in 
predator-prey dynamics: a non-parametric approach", Proceedings of 
the Royal Society of London B, 267, 1611-1620.
}
\examples{
data(veilleux)
## maybe str(veilleux) ; plot(veilleux) ...
}
\keyword{datasets}
