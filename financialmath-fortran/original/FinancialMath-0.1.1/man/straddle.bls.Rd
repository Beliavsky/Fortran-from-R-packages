\arguments{
  \item{S}{spot price at time 0}
  \item{K}{strike price of the call and put}
  \item{r}{continuously compounded yearly risk free rate}
  \item{t}{time of expiration (in years)}
  \item{sd}{standard deviation of the stock (volatility)}
  \item{position}{either buyer or seller of option ("long" or "short")}
  \item{plot}{tells whether or not to plot the payoff and profit}
}
\description{Gives a table and graphical representation of the payoff and profit of a long or short straddle for a range of future stock prices. Uses the Black Scholes equation for the call and put prices.}
\name{straddle.bls}
\alias{straddle.bls}
\title{Straddle Spread - Black Scholes}
\keyword{option}
\keyword{spread}
\keyword{put}
\keyword{call}
\value{A list of two components.
\item{Payoff}{A data frame of different payoffs and profits for given stock prices.}
\item{Premiums}{A matrix of the premiums for the call and put options, and the net cost.}
}
\examples{straddle.bls(S=100,K=110,r=.03,t=1,sd=.2,position="short")

straddle.bls(S=100,K=110,r=.03,t=1,sd=.2,position="long",plot=TRUE)}
\usage{straddle.bls(S,K,r,t,sd,position,plot=FALSE)}
\seealso{
  \code{\link{option.put}}

  \code{\link{option.call}}

  \code{\link{strangle.bls}}
}
\details{
Stock price at time t \eqn{=S_t}

Long Position:

For \eqn{S_t<=K}: payoff \eqn{=K-S_t}

For \eqn{S_t>K}: payoff \eqn{=S_t-K}

profit = payoff\eqn{-(price_{call}+price_{put})*e^{r*t}}

Short Position:

For \eqn{S_t<=K}: payoff \eqn{=S_t-K}

For \eqn{S_t>K}: payoff \eqn{=K-S_t}

profit = payoff\eqn{+(price_{call}+price_{put})*e^{r*t}}
}
