\name{bls.order1}
\alias{bls.order1}
\title{Black Scholes First-order Greeks}
\description{Gives the price and first order greeks for call and put options in the Black Scholes equation.}
\usage{bls.order1(S,K,r,t,sd,D=0)}
\author{Kameron Penn and Jack Schmidt}
\arguments{
  \item{S}{spot price at time 0}
  \item{K}{strike price}
  \item{r}{continuously compounded yearly risk free rate}
  \item{t}{time of expiration (in years)}
  \item{sd}{standard deviation of the stock (volatility)}
  \item{D}{continuous dividend yield}
}
\examples{
x <- bls.order1(S=100, K=110, r=.05, t=1, sd=.1, D=0)
ThetaPut <- x["Theta","Put"]
DeltaCall <- x[2,1]
}
\note{Cannot have any inputs as vectors.

t cannot be negative.

Either both or neither of S and K must be negative.}
\keyword{greeks}
\keyword{call}
\keyword{put}
\keyword{option}
\value{A matrix of the calculated greeks and prices for call and put options.}
\seealso{
  \code{\link{option.put}}

  \code{\link{option.call}}
}
