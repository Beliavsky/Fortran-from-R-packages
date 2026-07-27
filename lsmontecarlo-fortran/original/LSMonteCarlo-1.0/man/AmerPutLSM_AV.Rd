\name{AmerPutLSM_AV}
\alias{AmerPutLSM_AV}
\alias{print.AmerPutAV}
\alias{summary.AmerPutAV}

\title{
Pricing plain vanilla American put with Antithetic Variates }

\description{
The function calculates the price of a plain vanilla American put with Least Squares Monte Carlo method with Antithetic Variates (Glasserman, 2004). The regression model included in the algorithm is quadratic polynomial (Longstaff & Schwartz, 2000).
}
\usage{
AmerPutLSM_AV(Spot = 1, sigma = 0.2, n = 1000, m = 365, Strike = 1.1, r = 0.06, 
dr = 0, mT = 1)

\method{print}{AmerPutAV}(x, ...)
\method{summary}{AmerPutAV}(object, ...)

}
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
  \item{Strike}{
Strike price of the option.
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

\item{x}{
An object returned by the functions \code{AmerPutLSM_AV}. }
\item{object}{
An object returned by the function \code{AmerPutLSM_AV}. }
\item{...}{Not used. }
}
\value{
The function returns an object of the class AmerPutAV that is a list comprising the price calculated and the entry parameters. Class-specific \code{print} function gives the option type information and the price. The price as a single number can be derived using the \code{price} function. An overview of the entire object can be seen using the \code{summary} function.
}
\references{
Glasserman, P. 2004. Monte Carlo Methods in Financial Engineering. Springer.

Longstaff, F.A., and E.S. Schwartz. 2000. Valuing american option by simulation: A simple least-squared approach. The Review of Financial Studies. 14:113-147.
}
\author{
Mikhail A. Beketov
}
\seealso{
Functions: \code{\link{price}}, 
\code{\link{AmerPutLSM}}, 
\code{\link{AmerPutLSM_CV}}, 
\code{\link{AsianAmerPutLSM}}, and 
\code{\link{QuantoAmerPutLSM}}.
}
\examples{
AmerPutLSM_AV(n=500, m=50)
put<-AmerPutLSM_AV(Spot=14.2, Strike=16.5, n=200, m=50)
put
summary(put)
price(put)
put$price
}

\keyword{ Monte Carlo }
\keyword{ Option pricing }
\keyword{ American put }