\name{rv}
\alias{rv}
\title{
Residual variance from a locfit object.
}
\usage{
rv(fit)
}
\description{
  As part of the \code{\link{locfit}} fitting procedure, an estimate
of the residual variance is computed; the \code{rv} function extracts
the variance from the \code{"locfit"} object.
The estimate used is the residual sum of squares
(or residual deviance, for quasi-likelihood models),
divided by the residual degrees of freedom.

For likelihood (not quasi-likelihood) models, the estimate is 1.0.
}
\arguments{
\item{fit}{
\code{"locfit"} object.
}}
\value{
  Returns the residual variance estimate from the \code{"locfit"} object.
}
\examples{
data(ethanol)
fit <- locfit(NOx~E,data=ethanol)
rv(fit)
}
\seealso{
\link{locfit}, \link{rv<-}
}
\keyword{smooth}
% Converted by Sd2Rd version 0.2-a5.
