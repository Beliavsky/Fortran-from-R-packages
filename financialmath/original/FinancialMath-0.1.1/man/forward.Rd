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
\description{Gives a table and graphical representation of the payoff of a forward contract, and calculates the forward price for the contract.}
\name{forward}
\alias{forward}
\title{Forward Contract}
\keyword{forward}
\note{Leave an input variable as NA if it is not needed (ie. k=NA if div.structure="none").}
\value{A list of two components.
\item{Payoff}{A data frame of different payoffs for given stock prices.}
\item{Price}{The forward price of the contract.}
}
\examples{
forward(S=100,t=2,r=.03,position="short",div.structure="none")

forward(S=100,t=2,r=.03,position="long",div.structure="discrete",dividend=3,k=.02)

forward(S=100,t=1,r=.03,position="long",div.structure="continuous",D=.01)}
\usage{forward(S,t,r,position,div.structure="none",dividend=NA,df=1,D=NA,k=NA,plot=FALSE)}
\seealso{
  \code{\link{forward.prepaid}}
}
\details{
Stock price at time t \eqn{=S_t}

Long Position: payoff = \eqn{S_t} - forward price

Short Position: payoff = forward price - \eqn{S_t}

\bold{If div.structure = "none"}

forward price\eqn{=S*e^{r*t}}

\bold{If div.structure = "discrete"}

\eqn{eff.i=e^r-1}

\eqn{j=(1+eff.i)^{\frac{1}{df}}-1}

Number of dividends: \eqn{t^*=t*df}

if k = NA: forward price \eqn{=S*e^{r*t}-}dividend\eqn{*{s_{\left. {\overline {\, t^* \,}}\! \right |j}}}

if k != j: forward price \eqn{=S*e^{r*t}-}dividend\eqn{*\frac{1-(\frac{1+k}{1+j})^{t^*}}{j-k}*e^{r*t}}

if k = j: forward price \eqn{=S*e^{r*t}-}dividend\eqn{*\frac{t^*}{1+j}*e^{r*t}}

\bold{If div.structure = "continuous"}

forward price\eqn{=S*e^{(r-D)*t}}

}
