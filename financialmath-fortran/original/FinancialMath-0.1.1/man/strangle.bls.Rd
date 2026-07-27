\arguments{
  \item{S}{spot price at time 0}
  \item{K1}{strike price of the long put}
  \item{K2}{strike price of the long call}
  \item{r}{continuously compounded yearly risk free rate}
  \item{t}{time of expiration (in years)}
  \item{sd}{standard deviation of the stock (volatility)}
  \item{plot}{tells whether or not to plot the payoff and profit}
}
\name{strangle.bls}
\alias{strangle.bls}
\title{Strangle Spread - Black Scholes}
\description{Gives a table and graphical representation of the payoff and profit of a long strangle spread for a range of future stock prices. Uses the Black Scholes equation for the call prices.}
\usage{strangle.bls(S,K1,K2,r,t,sd,plot=FALSE)}
\author{Kameron Penn and Jack Schmidt}
\keyword{option}
\keyword{spread}
\keyword{call}
\keyword{put}
\note{K1 < S < K2 must be true.}
\value{A list of two components.
\item{Payoff}{A data frame of different payoffs and profits for given stock prices.}
\item{Premiums}{A matrix of the premiums for the call and put options, and the net cost.}
}
\examples{strangle.bls(S=105,K1=100,K2=110,r=.03,t=1,sd=.2)

strangle.bls(S=115,K1=50,K2=130,r=.03,t=1,sd=.2)}
\seealso{
  \code{\link{option.put}}

  \code{\link{option.call}}

  \code{\link{straddle.bls}}
}
\details{
Stock price at time t \eqn{=S_t}

For \eqn{S_t<=K1}: payoff \eqn{=K1-S_t}

For \eqn{K1<S_t<K2}: payoff \eqn{=0}

For \eqn{S_t>=K2}: payoff \eqn{=S_t-K2}

profit = payoff\eqn{-(price_{K1}+price_{K2})*e^{r*t}}
}
