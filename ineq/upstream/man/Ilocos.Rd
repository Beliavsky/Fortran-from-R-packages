\name{Ilocos}
\alias{Ilocos}
\title{Income Metadata from Ilocos, Philippines}
\description{Income metadata from surveys conducted by the
Philippines' National Statistics Office.}
\usage{data(Ilocos)}
\format{A data frame with 632 observations of 8 variables.
\describe{
\item{income}{total income of household,}
\item{sex}{sex of household head (\code{"male"} or \code{"female"}),}
\item{family.size}{family size (sometimes averaged over two semesters),}
\item{urbanity}{factor with levels \code{"rural"} and \code{"urban"},}
\item{province}{factor indicating the particular province,}
\item{AP.income}{total household income during the APIS,}
\item{AP.family.size}{family size during APIS,}
\item{AP.weight}{APIS survey weight for each household.}}}

\details{The data contains household income and metadata in one
of the sixteen regions of the Philippines called Ilocos.
The data comes from two of the NSO's surveys:
the 1997 Family and Income and Expenditure Survey
and the 1998 Annual Poverty Indicators Survey (APIS).

Since the APIS only has a six month reference
period, the original data were rescaled using an adjustment factor
from the quarterly GDP figures that can be obtained from the
major sectors.}
\source{National Statistics Office, Philippines:
\url{http://www.census.gov.ph/}, where also the whole data set may be
obtained.}
\keyword{datasets}
