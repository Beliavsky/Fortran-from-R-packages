\arguments{
  \item{S}{spot price at time 0}
  \item{K}{strike price}
  \item{r}{continuously compounded yearly risk free rate}
  \item{t}{time of expiration (in years)}
  \item{sd}{standard deviation of the stock (volatility)}
  \item{price}{specified call price if the Black Scholes pricing is not desired (leave as NA to use the Black Scholes pricing)}
  \item{position}{either buyer or seller of option ("long" or "short")}
  \item{plot}{tells whether or not to plot the payoff and profit}
}
\name{option.call}
\alias{option.call}
\title{Call Option}
\description{Gives a table and graphical representation of the payoff and profit of a long or short call option for a range of future stock prices.}
\usage{option.call(S,K,r,t,sd,price=NA,position,plot=FALSE)}
\author{Kameron Penn and Jack Schmidt}
\keyword{call}
\keyword{option}
\note{Finds the call price by using the Black Scholes equation by default.}
\value{A list of two components.
\item{Payoff}{A data frame of different payoffs and profits for given stock prices.}
\item{Premium}{The price for the call option.}
}
\examples{option.call(S=100,K=110,r=.03,t=1.5,sd=.2,price=NA,position="short")

option.call(S=100,K=100,r=.03,t=1,sd=.2,price=10,position="long")}
\seealso{
  \code{\link{option.put}}

  \code{\link{bls.order1}}
}
\details{
Stock price at time t \eqn{=S_t}

Long Position:

payoff = max\eqn{(0,S_t-K)}

profit = payoff - price\eqn{*e^{r*t}}

Short Position:

payoff = -max\eqn{(0,S_t-K)}

profit = payoff + price\eqn{*e^{r*t}}
}
