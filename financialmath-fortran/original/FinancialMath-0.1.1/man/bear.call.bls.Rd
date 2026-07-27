\name{bear.call.bls}
\alias{bear.call.bls}
\title{Bear Call Spread - Black Scholes}
\description{Gives a table and graphical representation of the payoff and profit of a bear call spread for a range of future stock prices. Uses the Black Scholes equation for the call prices.}
\usage{bear.call.bls(S,K1,K2,r,t,sd,plot=FALSE)}
\author{Kameron Penn and Jack Schmidt}
\arguments{
  \item{S}{spot price at time 0}
  \item{K1}{strike price of the short call}
  \item{K2}{strike price of the long call}
  \item{r}{yearly continuously compounded risk free rate}
  \item{t}{time of expiration (in years)}
  \item{sd}{standard deviation of the stock (volatility)}
  \item{plot}{tells whether or not to plot the payoff and profit}
}
\examples{
bear.call.bls(S=100,K1=70,K2=130,r=.03,t=1,sd=.2)}
\note{K1 must be less than S, and K2 must be greater than S.}
\keyword{spread}
\keyword{call}
\keyword{option}
\value{A list of two components.
\item{Payoff}{A data frame of different payoffs and profits for given stock prices.}
\item{Premiums}{A matrix of the premiums for the call options and the net cost.}
}
\seealso{
  \code{\link{bear.call}}

  \code{\link{bull.call.bls}}

  \code{\link{option.call}}
}
\details{
Stock price at time t \eqn{=S_t}

For \eqn{S_t<=K1}: payoff \eqn{=0}

For \eqn{K1<S_t<K2}: payoff \eqn{=K1-S_t}

For \eqn{S_t>=K2}: payoff \eqn{=K1-K2}

payoff = profit\eqn{+(price_{K1}-price_{K2})*e^{r*t}}
}
