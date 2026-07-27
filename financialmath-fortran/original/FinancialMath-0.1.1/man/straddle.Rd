\arguments{
  \item{S}{spot price at time 0}
  \item{K}{strike price of the call and put}
  \item{r}{continuously compounded yearly risk free rate}
  \item{t}{time of expiration (in years)}
  \item{price1}{price of the long call with strike price K}
  \item{price2}{price of the long put with strike price K}
  \item{position}{either buyer or seller of option ("long" or "short")}
  \item{plot}{tells whether or not to plot the payoff and profit}
}
\description{Gives a table and graphical representation of the payoff and profit of a long or short straddle for a range of future stock prices. }
\name{straddle}
\alias{straddle}
\title{Straddle Spread}
\keyword{option}
\keyword{spread}
\keyword{put}
\keyword{call}
\value{A list of two components.
\item{Payoff}{A data frame of different payoffs and profits for given stock prices.}
\item{Premiums}{A matrix of the premiums for the call and put options, and the net cost.}
}
\examples{straddle(S=100,K=110,r=.03,t=1,price1=15,price2=10,position="short")}
\usage{straddle(S,K,r,t,price1,price2,position,plot=FALSE)}
\seealso{
  \code{\link{straddle.bls}}

  \code{\link{option.put}}

  \code{\link{option.call}}

  \code{\link{strangle}}
}
\details{
Stock price at time t \eqn{=S_t}

Long Position:

For \eqn{S_t<=K}: payoff \eqn{=K-S_t}

For \eqn{S_t>K}: payoff \eqn{=S_t-K}

profit = payoff - (price1 + price2)\eqn{*e^{r*t}}

Short Position:

For \eqn{S_t<=K}: payoff \eqn{=S_t-K}

For \eqn{S_t>K}: payoff \eqn{=K-S_t}

profit = payoff + (price1 + price2)\eqn{*e^{r*t}}
}
