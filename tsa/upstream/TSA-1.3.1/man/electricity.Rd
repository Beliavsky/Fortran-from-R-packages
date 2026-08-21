\name{electricity}
\alias{electricity}
\docType{data}
\title{ Monthly US electricity production / time series}
\description{
Monthly U.S. electricity generation (in millions of kilowatt hours)
of all types: coal, natural gas, nuclear, petroleum, and wind,
01/1973 - 12/2005.
}
\usage{data(electricity)}
\format{
  The format is:
 'ts' int [1:396, 1] 160218 143539 148158 139589 147395 161244 173733 177365 156875 154197 ...
 - attr(*, "dimnames")=List of 2
  ..$ : NULL
  ..$ : chr "electricity"
 - attr(*, "tsp")= num [1:3] 1973 2006   12
}
\source{
Source: www.eia.doe.gov/emeu/mer/elect.html
}
\examples{
data(electricity)
## maybe str(electricity) ; plot(electricity) ...
}
\keyword{datasets}
