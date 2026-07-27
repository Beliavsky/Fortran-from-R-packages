\arguments{
  \item{S}{spot price at time 0}
  \item{K1}{strike price of the long put}
  \item{K2}{strike price of the short call}
  \item{r}{yearly continuously compounded risk free rate}
  \item{t}{time of expiration (in years)}
  \item{sd}{standard deviation of the stock (volatility)}
  \item{plot}{tells whether or not to plot the payoff and profit}
}
\description{Gives a table and graphical representation of the payoff and profit of a collar strategy for a range of future stock prices. Uses the Black Scholes equation for the call and put prices.}
\name{collar.bls}
\title{Collar Strategy - Black Scholes}
\alias{collar.bls}
\keyword{option}
\keyword{call}
\keyword{put}
\keyword{spread}
\value{A list of two components.
\item{Payoff}{A data frame of different payoffs and profits for given stock prices.}
\item{Premiums}{A matrix of the premiums for the call and put options and the net cost.}
}
\examples{collar.bls(S=100,K1=90,K2=110,r=.05,t=1,sd=.2)}
\usage{collar.bls(S,K1,K2,r,t,sd,plot=FALSE)}
\seealso{
  \code{\link{option.put}}

  \code{\link{option.call}}
}
\details{
Stock price at time t \eqn{=S_t}

For \eqn{S_t<=K1}: payoff \eqn{=K1-S_t}

For \eqn{K1<S_t<K2}: payoff \eqn{=0}

For \eqn{S_t>=K2}: payoff \eqn{=K2-S_t}

profit = payoff\eqn{+(price_{K2}-price_{K1})*e^{r*t}}
}
