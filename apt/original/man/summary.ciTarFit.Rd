\name{summary.ciTarFit}
\alias{summary.ciTarFit}
\title{Summary of Results from Threshold Cointegration Regression}
\description{This summarizes the main results from threshold cointegration regression.}
\usage{\method{summary}{ciTarFit}(object, digits=3, \dots)}
\arguments{
  \item{object}{an object of class 'ciTarFit'.}
  \item{digits}{number of digits for rounding outputs.}
  \item{\dots}{additional arguments to be passed.}
}
\details{
This wraps up the outputs from threshold cointegration regression in two data frames, one for diagnostic statistics and the other for coeffients.
}
\value{A list with two data frames. \code{dia} contains the main model specifications and hypothesis test results. \code{out} contains the regression results for both the long run (LR) and threshold cointegration (CI).
}
\author{Changyou Sun (\email{edwinsun258@gmail.com})}

\seealso{\code{\link{ciTarFit}}.}
\examples{
# see example at daVich
} 
\keyword{methods}
