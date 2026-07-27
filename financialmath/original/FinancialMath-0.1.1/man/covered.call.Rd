\arguments{
  \item{S}{spot price at time 0}
  \item{K}{strike price}
  \item{r}{continuously compounded yearly risk free rate}
  \item{t}{time of expiration (in years)}
  \item{sd}{standard deviation of the stock (volatility)}
  \item{price}{specified call price if the Black Scholes pricing is not desired (leave as NA to use the Black Scholes pricing)}
  \item{plot}{tells whether or not to plot the payoff and profit}
}
\description{Gives a table and graphical representation of the payoff and profit of a covered call strategy for a range of future stock prices.}
\alias{covered.call}
\name{covered.call}
\title{Covered Call}
\keyword{call}
\keyword{option}
\note{Finds the put price by using the Black Scholes equation by default.}
\value{A list of two components.
\item{Payoff}{A data frame of different payoffs and profits for given stock prices.}
\item{Premium}{The price of the call option.}
}
\examples{covered.call(S=100,K=110,r=.03,t=1,sd=.2,plot=TRUE)}
\usage{covered.call(S,K,r,t,sd,price=NA,plot=FALSE)}
\seealso{
  \code{\link{option.call}}

  \code{\link{covered.put}}
}
\details{
Stock price at time t \eqn{=S_t}

For \eqn{S_t<=K}: payoff \eqn{=S_t}

For \eqn{S_t>K}: payoff \eqn{=K}

profit = payoff + price\eqn{*e^{r*t}-S}
}
