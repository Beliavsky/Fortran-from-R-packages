\arguments{
  \item{S}{spot price at time 0}
  \item{K1}{strike price of the long put}
  \item{K2}{strike price of the short call}
  \item{r}{yearly continuously compounded risk free rate}
  \item{t}{time of expiration (in years)}
  \item{price1}{price of the long put with strike price K1}
  \item{price2}{price of the short call with strike price K2}
  \item{plot}{tells whether or not to plot the payoff and profit}
}
\description{Gives a table and graphical representation of the payoff and profit of a collar strategy for a range of future stock prices.}
\name{collar}
\title{Collar Strategy}
\alias{collar}
\keyword{option}
\keyword{call}
\keyword{put}
\keyword{spread}
\value{A list of two components.
\item{Payoff}{A data frame of different payoffs and profits for given stock prices.}
\item{Premiums}{A matrix of the premiums for the call and put options and the net cost.}
}
\examples{collar(S=100,K1=90,K2=110,r=.05,t=1,price1=5,price2=15,plot=TRUE)}
\usage{collar(S,K1,K2,r,t,price1,price2,plot=FALSE)}
\seealso{
  \code{\link{collar.bls}}

  \code{\link{option.put}}

  \code{\link{option.call}}
}
\details{
Stock price at time t \eqn{=S_t}

For \eqn{S_t<=K1}: payoff \eqn{=K1-S_t}

For \eqn{K1<S_t<K2}: payoff \eqn{=0}

For \eqn{S_t>=K2}: payoff \eqn{=K2-S_t}

profit = payoff + (price2 - price1)\eqn{*e^{r*t}}
}
