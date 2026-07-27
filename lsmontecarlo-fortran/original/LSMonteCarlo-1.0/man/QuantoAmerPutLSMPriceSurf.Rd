\name{QuantoAmerPutLSMPriceSurf}
\alias{QuantoAmerPutLSMPriceSurf}

\title{
Deriving a table of Quanto American put prices at different volatilities and strikes
}
\description{
The function calculates the prices at different volatilities and strikes using the \code{QuantoAmerPutLSM} function.
}
\usage{
QuantoAmerPutLSMPriceSurf(Spot = 1, vols = (seq(0.1, 2, 0.1)), n = 1000, m = 365, 
strikes = (seq(0.5, 2.5, 0.1)), r = 0.06, dr = 0, mT = 1, Spot2 = 1, sigma2 = 0.2, 
r2 = 0, dr2 = 0, rho = 0)
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
Interest rate of the numeraire currency (e.g. USD).
}
  \item{dr}{
Dividend rate of the underlying asset.
}
  \item{mT}{
Maturity time (years).
}
  \item{Spot2}{
Spot price of the 3rd asset (e.g. EUR/USD).
}
  \item{sigma2}{
Volatility of the 3rd asset.
}
  \item{r2}{
Interest rate of the 3rd asset.
}
  \item{dr2}{
Dividend rate of the 3rd asset.
}
  \item{rho}{
Correlation coefficient between the prices. 
}
}

\value{
The function returns an object of the class PriceSurface that is a matrix of prices at different volatilities and strikes. Class-specific \code{summary} function gives the sequences of volatilities and strikes used, as well as maximum, minimum, and average prices. Class-specific \code{plot} function constructs a 3-D plot of the price surface.
}

\author{
Mikhail A. Beketov
}
\note{
The function \code{rmvnorm} included in the pricing algorithm is a part of the mnormt package. Please, load that package before the use of the \code{QuantoAmerPutLSMPriceSurf} function. Using the function \code{plot} with default pallet requires package fBasics, if you do not want to load this package, you can set \code{color=NULL} or other palette).
}


\seealso{
Functions: \code{\link{QuantoAmerPutLSM}},
\code{\link{summary.PriceSurface}}, 
\code{\link{plot.PriceSurface}}, 
\code{\link{AmerPutLSMPriceSurf}}, and 
\code{\link{AsianAmerPutLSMPriceSurf}}.
}
\examples{
surface<-QuantoAmerPutLSMPriceSurf(vols = (seq(0.1, 1.7, 0.2)), n=100, m=5, 
strikes = (seq(0.7, 1.7, 0.2)))
summary(surface)
plot(surface, color = divPalette(150, "RdBu"))
}

\keyword{ Monte Carlo }
\keyword{ Option pricing }
\keyword{ American put }
