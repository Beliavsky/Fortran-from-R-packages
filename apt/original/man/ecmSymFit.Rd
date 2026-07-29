\name{ecmSymFit}
\alias{ecmSymFit}
\title{Fitting symmetric Error Correction Model}
\description{Estimate a symmetric error correction model (ECM) for two time series.}
\usage{
ecmSymFit(y, x, lag = 1)
}
\arguments{
  \item{y}{dependent or left-side variable for the long-run regression.}
  \item{x}{independent or right-side variable for the long-run regression.}
  \item{lag}{number of lags for variables on the right side.}
}
\details{
The package \code{apt} focuses on price transmission between two series. This function estimates a standard error correction model for two time series. While it can be extended for more than two series, it is beyond the objective of the package now. 
}
\value{
Return a list object of class "ecm" and "ecmSymFit" with the following components: 
  \item{y}{dependend variable} 
  \item{x}{independent variable}
  \item{lag}{number of lags}    
  \item{data}{all the data combined for the ECM}
  \item{IndVar}{data frame of the right-hand variables used in the ECM}
  \item{name.x}{name of the independent variable}
  \item{name.y}{name of the dependent variable}  
  \item{ecm.y}{ECM regression for the dependent variable }
  \item{ecm.x}{ECM regression for the independent variable}
}
\references{
Enders, W. 2004. Applied Econometric Time Series. John Wiley & Sons, Inc., New York. 480 P.
}
\author{Changyou Sun (\email{edwinsun258@gmail.com})}
\seealso{
\code{\link{print.ecm}}; \code{\link{summary.ecm}}; \code{\link{ecmDiag}}; and \code{\link{ecmAsyFit}}.
}
\examples{
# see example at daVich
}
\keyword{regression}