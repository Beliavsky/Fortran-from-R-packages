\arguments{
  \item{prices}{vector of variable prices}
  \item{rates}{vector of variable rates}
  \item{type}{rates defined as either "spot_rate" or "zcb_price"}
}
\name{swap.commodity}
\title{Commodity Swap}
\alias{swap.commodity}
\description{Solves for the fixed swap price, given the variable prices and interest rates (either as spot rates or zero coupon bond prices).}
\usage{swap.commodity(prices, rates, type="spot_rate")}
\author{Kameron Penn and Jack Schmidt}
\keyword{swap}
\note{Length of the price vector and rate vector must be of the same length.
}
\value{The fixed swap price.}
\examples{swap.commodity(prices=c(103,106,108), rates=c(.04,.05,.06))

swap.commodity(prices=c(103,106,108), rates=c(.9615,.907,.8396),type="zcb_price")

swap.commodity(prices=c(105,105,105), rates=c(.85,.89,.80),type="zcb_price")
}
\seealso{
  \code{\link{swap.rate}}
}
\details{
For spot rates: \eqn{\sum_{k=1}^n\frac{prices_k}{(1+rates_k)^k}=\sum_{k=1}^n\frac{X}{(1+rates_k)^k}}

For zero coupon bond prices: \eqn{\sum_{k=1}^nprices_k*rates_k=\sum_{k=1}^nX*rates_k}

Where \eqn{X=} fixed swap price.
}
