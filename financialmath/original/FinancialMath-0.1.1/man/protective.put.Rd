\arguments{
  \item{S}{spot price at time 0}
  \item{K}{strike price}
  \item{r}{continuously compounded yearly risk free rate}
  \item{t}{time of expiration (in years)}
  \item{sd}{standard deviation of the stock (volatility)}
  \item{price}{specified put price if the Black Scholes pricing is not desired (leave as NA to use the Black Scholes pricing)}
  \item{plot}{tells whether or not to plot the payoff and profit}
}
\name{protective.put}
\alias{protective.put}
\title{Protective Put}
\description{Gives a table and graphical representation of the payoff and profit of a protective put strategy for a range of future stock prices.}
\usage{protective.put(S,K,r,t,sd,price=NA,plot=FALSE)}
\author{Kameron Penn and Jack Schmidt}
\keyword{put}
\keyword{option}
\note{Finds the put price by using the Black Scholes equation by default.}
\value{A list of two components.
\item{Payoff}{A data frame of different payoffs and profits for given stock prices.}
\item{Premium}{The price of the put option.}
}
\examples{protective.put(S=100,K=100,r=.03,t=1,sd=.2)

protective.put(S=100,K=90,r=.01,t=.5,sd=.1)}
\seealso{
  \code{\link{option.put}}
}
\details{
Stock price at time t \eqn{=S_t}

For \eqn{S_t<=K}: payoff \eqn{=K-S}

For \eqn{S_t>K}: payoff \eqn{=S_t-S}

profit = payoff - price\eqn{*e^{r*t}}
}
