\name{fastGBM}
\alias{fastGBM}

\title{
Generating Geometric Brownian motion
}
\description{
Quick Generating Geometric Brownian motion avoiding unnecessary loops using the \code{cumsum} function. Technical function implemented in the pricing functions of the package.
}
\usage{
fastGBM(Spot = 1, sigma = 0.2, n = 1000, m = 365, r = 0.06, dr = 0, mT = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Spot}{
Spot price of the underlying asset (e.g. stock).
}
  \item{sigma}{
Volatility of the underlying asset.
}
  \item{n}{
Number of paths simulated.
}
  \item{m}{
Number of time steps in the simulation.
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
}

\value{
Table with paths generated (each row is a path, class "matrix")  
}
\author{
Mikhail A. Beketov
}
\seealso{
Functions: \code{\link{AmerPutLSM}},  
\code{\link{AmerPutLSM_AV}}, 
\code{\link{AmerPutLSM_CV}}, 
\code{\link{AsianAmerPutLSM}}, 
\code{\link{QuantoAmerPutLSM}}, and 
\code{\link{QuantoAmerPutLSM_AV}}.
}
\examples{
fastGBM(n=10, m=5)
matplot(t(fastGBM(n=100, m=100)), type="l") # matrix transpose by "t()"
}

\keyword{ Monte Carlo }
\keyword{ Option pricing }
\keyword{ American put }
