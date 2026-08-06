\name{salmon}
\alias{salmon}
\docType{data}
\title{  Salmon data }
\description{
  The salmon data contains two measurements of the growth rings on 
  the scale of Alaskan and Canadian salmon as well as the gender 
  of the fishes. There are 50 Alaskan-born and 50 Canadian-born salmon, and this information
  is coded in the variable \code{Origin}.
}
\usage{data(salmon)}
\format{
  A data frame with 100 observations on the following 4 variables.
  \describe{
    \item{\code{Gender}}{female=1 and male=2}
    \item{\code{Freshwater}}{diameter of rings for the first-year freshwater growth (hundrets of an inch)}
    \item{\code{Marine}}{diameter of rings for the first-year marine growth (hundrets of an inch)}
    \item{\code{Origin}}{Origin of the fish: a factor with levels \code{Alaskan} \code{Canadian}}
  }
}

%\details{}

\source{
Johnson, R.A. and Wichern, D. W. \emph{Applied Multivariate
Statistical Analysis} (Prentice Hall, International
Editions, 2002, fifth edition)}

%\references{}
\examples{
data(salmon)
}
\keyword{datasets}
