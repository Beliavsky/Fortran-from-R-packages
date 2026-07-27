\arguments{
  \item{S}{spot price at time 0}
  \item{K1}{strike price of the long put}
  \item{K2}{strike price of the long call}
  \item{r}{continuously compounded yearly risk free rate}
  \item{t}{time of expiration (in years)}
  \item{price1}{price of the long put with strike price K1}
  \item{price2}{price of the long call with strike price K2}
  \item{plot}{tells whether or not to plot the payoff and profit}
}
\name{strangle}
\alias{strangle}
\title{Strangle Spread}
\description{Gives a table and graphical representation of the payoff and profit of a long strangle spread for a range of future stock prices.}
\usage{strangle(S,K1,K2,r,t,price1,price2,plot=FALSE)}
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
\examples{strangle(S=105,K1=100,K2=110,r=.03,t=1,price1=10,price2=15,plot=TRUE)}
\seealso{
  \code{\link{strangle.bls}}

  \code{\link{option.put}}

  \code{\link{option.call}}

  \code{\link{straddle}}
}
\details{
Stock price at time t \eqn{=S_t}

For \eqn{S_t<=K1}: payoff \eqn{=K1-S_t}

For \eqn{K1<S_t<K2}: payoff \eqn{=0}

For \eqn{S_t>=K2}: payoff \eqn{=S_t-K2}

profit = payoff - (price1 + price2)\eqn{*e^{r*t}}
}
