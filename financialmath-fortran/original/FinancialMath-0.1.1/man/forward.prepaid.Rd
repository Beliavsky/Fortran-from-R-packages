\arguments{
  \item{S}{spot price at time 0}
  \item{t}{time of expiration (in years)}
  \item{r}{continuously compounded yearly risk free rate}
  \item{position}{either buyer or seller of the contract ("long" or "short")}
  \item{div.structure}{the structure of the dividends for the underlying ("none", "continuous", or "discrete")}
  \item{dividend}{amount of each dividend, or amount of first dividend if k is not NA}
  \item{df}{dividend frequency- number of dividends per year}
  \item{D}{continuous dividend yield}
  \item{k}{dividend growth rate per df}
  \item{plot}{tells whether or not to plot the payoff}
}
\description{Gives a table and graphical representation of the payoff of a prepaid forward contract, and calculates the prepaid forward price for the contract.}
\name{forward.prepaid}
\alias{forward.prepaid}
\title{Prepaid Forward Contract}
\keyword{forward}
\note{Leave an input variable as NA if it is not needed (ie. k=NA if div.structure="none").}
\value{A list of two components.
\item{Payoff}{A data frame of different payoffs for given stock prices.}
\item{Price}{The prepaid forward price of the contract.}
}
\examples{
forward.prepaid(S=100,t=2,r=.04,position="short",div.structure="none")

forward.prepaid(S=100,t=2,r=.03,position="long",div.structure="discrete",
dividend=3,k=.02,df=2)

forward.prepaid(S=100,t=1,r=.05,position="long",div.structure="continuous",D=.06)}
\usage{forward.prepaid(S,t,r,position,div.structure="none",dividend=NA,df=1,D=NA,
k=NA,plot=FALSE)}
\seealso{
  \code{\link{forward}}
}
\details{
Stock price at time t \eqn{=S_t}

Long Position: payoff = \eqn{S_t} - prepaid forward price

Short Position: payoff = prepaid forward price - \eqn{S_t}

\bold{If div.structure = "none"}

forward price\eqn{=S}

\bold{If div.structure = "discrete"}

\eqn{eff.i=e^r-1}

\eqn{j=(1+eff.i)^{\frac{1}{df}}-1}

Number of dividends: \eqn{t^*=t*df}

if k = NA: prepaid forward price \eqn{=S-}dividend\eqn{*{a_{\left. {\overline {\, t^* \,}}\! \right |j}}}

if k != j: prepaid forward price \eqn{=S-}dividend\eqn{*\frac{1-(\frac{1+k}{1+j})^{t^*}}{j-k}}

if k = j: prepaid forward price \eqn{=S-}dividend\eqn{*\frac{t^*}{1+j}}

\bold{If div.structure = "continuous"}

prepaid forward price\eqn{=S*e^{-D*t}}

}
