\arguments{
  \item{f}{face value}
  \item{r}{coupon rate convertible cf times per year}
  \item{c}{redemption value}
  \item{n}{the number of coupons/periods for the bond}
  \item{i}{nominal interest rate convertible ic times per year}
  \item{ic}{interest conversion frequency per year}
  \item{cf}{coupon frequency- number of coupons per year}
  \item{t}{specified period for which the price and write up/down amount is solved for, if not NA}
  \item{plot}{tells whether or not to plot the convexity}
}
\description{Solves for the price, premium/discount, and Durations and Convexities (in terms of periods). At a specified period (t), it solves for the full and clean prices, and the write up/down amount. Also has the option to plot the convexity of the bond.}
\name{bond}
\alias{bond}
\title{Bond Analysis}
\keyword{bond}
\keyword{analysis}
\note{t must be less than n.

To make the duration in terms of years, divide it by cf.

To make the convexity in terms of years, divide it by \eqn{cf^2}.}
\value{A matrix of all of the bond details and calculated variables.}
\examples{
bond(f=100,r=.04,c=100,n=20,i=.04,ic=1,cf=1,t=1)

bond(f=100,r=.05,c=110,n=10,i=.06,ic=1,cf=2,t=5)}
\usage{bond(f,r,c,n,i,ic=1,cf=1,t=NA,plot=FALSE)}
\details{
Effective Rate of Interest: \eqn{eff.i=(1+\frac{i}{ic})^{ic}-1}

\eqn{j=(1+eff.i)^{\frac{1}{cf}}-1}

coupon \eqn{=\frac{f*r}{cf}} (per period)

price = coupon\eqn{*{a_{\left. {\overline {\, n \,}}\! \right |j}}+c*(1+j)^{-n}}

\eqn{MAC D=\frac{\sum_{k=1}^n k*(1+j)^{-k}*coupon+n*(1+j)^{-n}*c}{price}}

\eqn{MOD D=\frac{\sum_{k=1}^n k*(1+j)^{-(k+1)}*coupon+n*(1+j)^{-(n+1)}*c}{price}}

\eqn{MAC C=\frac{\sum_{k=1}^n k^2*(1+j)^{-k}*coupon+n^2*(1+j)^{-n}*c}{price}}

\eqn{MOD C=\frac{\sum_{k=1}^n k*(k+1)*(1+j)^{-(k+2)}*coupon+n*(n+1)*(1+j)^{-(n+2)}*c}{price}}

\bold{Price (for period t):}

If t is an integer: price =coupon\eqn{*{a_{\left. {\overline {\, n-t \,}}\! \right |j}}+c*(1+j)^{-(n-t)}}

If t is not an integer then \eqn{t=t^*+k} where \eqn{t^*} is an integer and \eqn{0<k<1}:

full price \eqn{=(} coupon\eqn{*{a_{\left. {\overline {\, n-t^* \,}}\! \right |j}}+c*(1+j)^{-(n-t^*)})*(1+j)^k}

clean price = full price\eqn{-k*}coupon

\bold{If price > c :}

premium = price\eqn{-c}

Write-down amount (for period t) \eqn{=(}coupon\eqn{-c*j)*(1+j)^{-(n-t+1)}}

\bold{If price < c :}

discount \eqn{=c-}price

Write-up amount (for period t) \eqn{=(c*j-}coupon\eqn{)*(1+j)^{-(n-t+1)}}

}
