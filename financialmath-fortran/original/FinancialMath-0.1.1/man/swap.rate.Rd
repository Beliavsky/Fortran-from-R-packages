\arguments{
  \item{rates}{vector of variable rates}
  \item{type}{rates as either "spot_rate" or "zcb_price"}
}
\description{Solves for the fixed interest rate given the variable interest rates (either as spot rates or zero coupon bond prices).}
\name{swap.rate}
\title{Interest Rate Swap}
\keyword{swap}
\keyword{interest}
\alias{swap.rate}
\value{The fixed interest rate swap.}
\examples{swap.rate(rates=c(.04, .05, .06), type = "spot_rate")

swap.rate(rates=c(.93,.95,.98,.90), type = "zcb_price")}
\usage{swap.rate(rates, type="spot_rate")}
\seealso{
  \code{\link{swap.commodity}}
}
\details{
For spot rates: \eqn{1=\sum_{k=1}^n[\frac{R}{(1+rates_k)^k}]+\frac{1}{(1+rates_n)^n}}

For zero coupon bond prices: \eqn{1=\sum_{k=1}^n(R*rates_k)+rates_n}

Where \eqn{R=} fixed swap rate.
}
