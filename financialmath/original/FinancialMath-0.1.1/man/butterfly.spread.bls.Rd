\arguments{
  \item{S}{spot price at time 0}
  \item{K1}{strike price of the first long call}
  \item{K2}{strike price of the two short calls}
  \item{K3}{strike price of the second long call}
  \item{r}{continuously compounded yearly risk free rate}
  \item{t}{time of expiration (in years)}
  \item{sd}{standard deviation of the stock (volatility)}
  \item{plot}{tells whether or not to plot the payoff and profit}
}
\description{Gives a table and graphical representation of the payoff and profit of a long butterfly spread for a range of future stock prices. Uses the Black Scholes equation for the call prices.}
\name{butterfly.spread.bls}
\title{Butterfly Spread - Black Scholes}
\alias{butterfly.spread.bls}
\keyword{spread}
\keyword{option}
\keyword{call}
\note{K2 must be equal to S.

K3 and K1 must both be equidistant to K2 and S.

K1 < K2 < K3 must be true.}
\value{A list of two components.
\item{Payoff}{A data frame of different payoffs and profits for given stock prices.}
\item{Premiums}{A matrix of the premiums for the call options and the net cost.}
}
\examples{butterfly.spread.bls(S=100,K1=75,K2=100,K3=125,r=.03,t=1,sd=.2)}
\usage{butterfly.spread.bls(S,K1,K2=S,K3,r,t,sd,plot=FALSE)}
\seealso{
  \code{\link{butterfly.spread}}

  \code{\link{option.call}}
}
\details{
Stock price at time t \eqn{=S_t}

For \eqn{S_t<=K1}: payoff \eqn{=0}

For \eqn{K1<S_t<=K2}: payoff \eqn{=S_t-K1}

For \eqn{K2<S_t<K3}: payoff \eqn{=2*K2-K1-S_t}

For \eqn{S_t>=K3}: payoff \eqn{=0}

profit = payoff\eqn{+(2*price_{K2}-price_{K1}-price_{K3})*e^{r*t}}
}
