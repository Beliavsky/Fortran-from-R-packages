\name{dat}
\alias{dat}
\title{
Locfit - data evaluation structure.
}
\usage{
dat(cv=FALSE)
}
\description{
  \code{dat} is used to specify evaluation on the given data points 
  for \code{\link{locfit.raw}()}. 
}
\arguments{
\item{cv}{Whether cross-validation should be done.}
}
%\examples{
%data(ethanol, package="locfit")
%plot.eval(locfit(NOx~E+C, data=ethanol, scale=0, ev=dat()))
%}
\keyword{smooth}
