\name{AsianCall_AppLord}
\alias{AsianCall_AppLord}
\title{
Asian Options - Approximation
}
\description{
The price of an arithmetic average Asian option is computed using
the approximation method of Lord. 
}
\usage{
AsianCall_AppLord(T=1, d=12, K=100, r=0.05, sigma=0.1, S0=100, all=TRUE)
}
\arguments{
  \item{T}{
\code{T} time to maturity (in years)
}
  \item{d}{
\code{d} number of controll points
}
  \item{K}{
\code{K} strike price
}
  \item{r}{
\code{r} risk free interest rate
}
  \item{sigma}{
\code{sigma} volatility (yearly)
}
  \item{S0}{
\code{S0} starting stockprice 
}
  \item{all}{
       \code{all} \code{TRUE} means that the full Asian Call option price is approximated
}
}
\details{
\code{AsianCall_AppLord()} uses a sophisticated approximation of Lord (2006).
}
\value{
  returns the approximate price.
}
\references{
Lord, R., Partially Exact and Bounded Approximations for Arithmetic Asian
Options, Journal of Computational Finance, Vol. 10, No. 2, pp. 1-52, 2006
}
\author{
Kemal Dingec, Wolfgang Hormann
}

\seealso{
\code{\link{OptionPricing-package}}
}
\examples{
AsianCall_AppLord(T = 1, d = 12, K = 100, r = 0.05, sigma = 0.25, S0 = 100, all = TRUE)
}