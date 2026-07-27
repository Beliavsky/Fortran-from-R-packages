\name{QuantoAmerPutLSM}
\alias{QuantoAmerPutLSM}
\alias{print.QuantoAmerPut}
\alias{summary.QuantoAmerPut}

\title{
Calculating the price of Quanto American put
}
\description{
The function calculates the price of Quanto American put with Least Squares Monte Carlo method. The Quanto option is cash-settled option, whose pay-off is converted into a third currency/asset at exercise at a pre-specified rate/price (Wystup, 2011), and can also be considered as a usual option but settled in a "wrong" asset (Vecer, 2011). The regression model included in the algorithm is quadratic polynomial (Longstaff & Schwartz, 2000).
}
\usage{
QuantoAmerPutLSM(Spot = 1, sigma = 0.2, n = 1000, m = 365, Strike = 1.1, r = 0.06, 
dr = 0, mT = 1, Spot2 = 1, sigma2 = 0.2, r2 = 0, dr2 = 0, rho = 0)

\method{print}{QuantoAmerPut}(x, ...)
\method{summary}{QuantoAmerPut}(object, ...)


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

\item{x}{
An object returned by the functions \code{QuantoAmerPutLSM}. }
\item{object}{
An object returned by the function \code{QuantoAmerPutLSM}. }
\item{...}{Not used. }

}

\value{
The function returns an object of the class QuantoAmerPut that is a list comprising the price calculated, option type, and the entry parameters. Class-specific \code{print} function gives the option type information and the price. The price as a single number can be derived using the \code{price} function. An overview of the entire object can be seen using the \code{summary} function.
}
\references{
Longstaff, F.A., and E.S. Schwartz. 2000. Valuing american option by simulation: A simple least-squared approach. The Review of Financial Studies. 14:113-147.

Vecer, J. 2011. Stochastic Finance: A Numeraire Approach. CRC Press.

Wystup, U. 2011. Quanto Options. MathFinance AG.
}
\author{
Mikhail A. Beketov
}
\note{
The function \code{rmvnorm} included in the pricing algorithm is a part of the mnormt package. Please, load that package before the use of the \code{QuantoAmerPutLSM} function.
}

\seealso{
Functions: \code{\link{price}}, 
\code{\link{QuantoAmerPutLSM_AV}}, 
\code{\link{AmerPutLSM}}, 
\code{\link{AsianAmerPutLSM}}, and 
\code{\link{AmerPutLSM_AV}}.
}
\examples{
QuantoAmerPutLSM(n=200, m=50)
put<-QuantoAmerPutLSM(Spot=14.2, Strike=16.5, n=200, m=50)
put
summary(put)
price(put)
}

\keyword{ Monte Carlo }
\keyword{ Option pricing }
\keyword{ American put }
