\name{SharpeRatio.annualized}
\alias{SharpeRatio.annualized}
\encoding{latin1}
\title{calculate annualized Sharpe Ratio}
\usage{
SharpeRatio.annualized(R, Rf = 0, alpha=0.05,scale = NA, geometric = TRUE, FUN = "StdDev")
}
\arguments{
\item{R}{an xts, vector, matrix, data frame, timeSeries or zoo object of
asset returns}

\item{Rf}{risk free rate, in same period as your returns}

\item{alpha}{Tail probability for VaR or ES, default alpha=.05}

\item{scale}{number of periods in a year (daily scale = 252, monthly scale =
12, quarterly scale = 4)}

\item{geometric}{utilize geometric chaining (TRUE) or simple/arithmetic chaining (FALSE) to aggregate returns,default TRUE}

\item{FUN}{one of "StdDev" or "VaR" or "ES" to use as the denominator, default="StdDev"}

}
\description{
The Sharpe Ratio is a risk-adjusted measure of return that uses standard
deviation to represent risk.
}
\details{
The Sharpe ratio is simply the return per unit of risk (represented by
variance).  The higher the Sharpe ratio, the better the combined performance
of "risk" and return.

This function annualizes the number based on the scale parameter.

\deqn{\frac{\sqrt[n]{prod(1+R_{a})^{scale}}-1}{\sqrt{scale}\cdot\sqrt{\sigma}}}

Using an annualized Sharpe Ratio is useful for comparison of multiple return
streams.  The annualized Sharpe ratio is computed by dividing the annualized
mean monthly excess return by the annualized monthly standard deviation of
excess return.

William Sharpe now recommends Information Ratio preferentially to the
original Sharpe Ratio.
}
\examples{

  data(assetReturns)
	R=assetReturns[, -29]
  SharpeRatio.annualized(R)

}
\references{
Sharpe, W.F. The Sharpe Ratio,\emph{Journal of Portfolio
Management},Fall 1994, 49-58.\cr
See also package \code{PerformanceAnalytics}.
}
\seealso{
\code{\link{SharpeRatio}} \cr \code{\link{InformationRatio}} \cr
\code{\link{TrackingError}} \cr \code{\link{ActivePremium}} \cr
\code{\link{SortinoRatio}}
}


\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

