\name{EuPutBS}
\alias{EuPutBS}
\alias{EuCallBS}
\title{
Black & Scholes solution for European put and call
}
\description{
Pricing plain vanilla American put and call options using Black & Scholes solution.
}
\usage{
EuPutBS(Spot, sigma, Strike, r, dr, mT)

EuCallBS(Spot, sigma, Strike, r, dr, mT)
}

\arguments{
  \item{Spot}{
Spot price of the underlying asset (e.g. stock).
}
  \item{sigma}{
Volatility of the underlying asset.
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
}

\value{
The function returns the price as a single number (class "numeric"). 
}

\seealso{
\code{\link{AmerPutLSM_CV}}
}
\examples{
EuPutBS(1, 0.2, 1, 0.06, 0, 1)
EuCallBS(1, 0.2, 1, 0.06, 0, 1)
}

\keyword{ Option pricing }
\keyword{ American put }
\keyword{ American call }
