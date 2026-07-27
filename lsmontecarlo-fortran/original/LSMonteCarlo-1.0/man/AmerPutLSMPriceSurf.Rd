\name{AmerPutLSMPriceSurf}
\alias{AmerPutLSMPriceSurf}
\alias{summary.PriceSurface}
\alias{plot.PriceSurface}

\title{
Deriving a table of American put prices at different volatilities and strikes
}
\description{
The function calculates the prices at different volatilities and strikes using the \code{AmerPutLSM} function.
}
\usage{
AmerPutLSMPriceSurf(Spot = 1, vols = (seq(0.1, 2, 0.1)), n = 1000, m = 365, 
strikes = (seq(0.5, 2.5, 0.1)), r = 0.06, dr = 0, mT = 1)

\method{summary}{PriceSurface}(object, ...)
\method{plot}{PriceSurface}(x, color = divPalette(800, "RdBu"), ...)
}

\arguments{
  \item{Spot}{
Spot price of the underlying asset (e.g. stock).
}
  \item{vols}{
Sequence of volatilities.
}
  \item{n}{
Number of paths simulated.
}
  \item{m}{
Number of time steps in the simulation.
}
  \item{strikes}{
Sequence of strikes.
}
  \item{r}{
Interest rate of the numeraire currency (e.g. EUR).
}
  \item{dr}{
Dividend rate of the underlying asset.
}
  \item{mT}{
Maturity time (years).
}

  \item{object}{
Object of the class PriceSurface that is a matrix of prices at different volatilities and strikes. }
  \item{x}{
Object of the class PriceSurface that is a matrix of prices at different volatilities and strikes. }
  \item{color}{
Color palette (the default pallet requires package fBasics, if you do not want to load this package, you can set \code{color=NULL} or other palette).
}
  \item{...}{Not used. }
}

\value{
The function returns an object of the class PriceSurface that is a matrix of prices at different volatilities and strikes. Class-specific \code{summary} function gives the sequences of volatilities and strikes used, as well as maximum, minimum, and average prices. Class-specific \code{plot} function constructs a 3-D plot of the price surface.
}

\author{
Mikhail A. Beketov
}

\seealso{
Functions: \code{\link{AmerPutLSM}},  
\code{\link{AsianAmerPutLSMPriceSurf}}, and 
\code{\link{QuantoAmerPutLSMPriceSurf}}.
}
\examples{
surface<-AmerPutLSMPriceSurf(vols = (seq(0.1, 1.5, 0.2)), n=200, m=10, 
strikes = (seq(0.5, 1.9, 0.2)))
summary(surface)
plot(surface, color = divPalette(150, "RdBu"))
}

\keyword{ Monte Carlo }
\keyword{ Option pricing }
\keyword{ American put }
