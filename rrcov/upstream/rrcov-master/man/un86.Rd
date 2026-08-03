\name{un86}
\alias{un86}
\docType{data}
\title{ United Nations Data - 1986}
\description{
  This data set consists of seven socioeconomic variables observed for 73 countries.
}
\usage{data(un86)}
\format{
  A data frame with 73 observations on the following 7 variables.
  \describe{
    \item{\code{POP}}{Total population in milions}
    \item{\code{MOR}}{Number of infant deaths per thousand births}
    \item{\code{CAR}}{Number of motorized vehicles per hundred inhabitants}
    \item{\code{DR}}{Number of medical doctors per thousand inhabitants}
    \item{\code{GNP}}{Gross national product per inhabitant in thousands of US dollars}
    \item{\code{DEN}}{Density in inhabitants per square kilometer}
    \item{\code{TB}}{Trade balance, defined as total exports/(total exports + total imports)}
  }
}
\details{
  The data set is from World Statistics in Brief, Number 10, a 1986 UN publication. It was used in Daigle et al. (1992) to illustrate a robust biplot method.
}
\source{
    World Statistics in Brief, Number 10, a 1986 United Nations publication 
    
    Daigle, G. and Rivest, L. (1992) A robust biplot, The canadian Journal of Statistics, 
    20, pp 241--255
}
\examples{
data(un86)
pairs(un86)
}
\keyword{datasets}
