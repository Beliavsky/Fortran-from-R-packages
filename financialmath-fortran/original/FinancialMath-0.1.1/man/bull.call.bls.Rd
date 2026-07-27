\arguments{
  \item{S}{spot price at time 0}
  \item{K1}{strike price of the long call}
  \item{K2}{strike price of the short call}
  \item{r}{yearly continuously compounded risk free rate}
  \item{t}{time of expiration (in years)}
  \item{sd}{standard deviation of the stock (volatility)}
  \item{plot}{tells whether or not to plot the payoff and profit}
}
\description{Gives a table and graphical representation of the payoff and profit of a bull call spread for a range of future stock prices. Uses the Black Scholes equation for the call prices.}
\name{bull.call.bls}
\title{Bull Call Spread - Black Scholes}
\alias{bull.call.bls}
\keyword{spread}
\keyword{call}
\keyword{option}
\note{K1 must be less than S, and K2 must be greater than S.}
\value{A list of two components.
\item{Payoff}{A data frame of different payoffs and profits for given stock prices.}
\item{Premiums}{A matrix of the premiums for the call options and the net cost.}
}
\examples{
bull.call.bls(S=115,K1=100,K2=145,r=.03,t=1,sd=.2)}
\usage{bull.call.bls(S,K1,K2,r,t,sd,plot=FALSE)}
\seealso{
  \code{\link{bear.call}}

  \code{\link{option.call}}
}
\details{
Stock price at time t \eqn{=S_t}

For \eqn{S_t<=K1}: payoff \eqn{=0}

For \eqn{K1<S_t<K2}: payoff \eqn{=S_t-K1}

For \eqn{S_t>=K2}: payoff \eqn{=K2-K1}

profit = payoff\eqn{+(price_{K2}-price_{K1})*e^{r*t}}
}
