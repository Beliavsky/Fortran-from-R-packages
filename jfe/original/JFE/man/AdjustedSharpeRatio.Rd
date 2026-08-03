\name{AdjustedSharpeRatio}
\alias{AdjustedSharpeRatio}
\encoding{latin1}

\title{Adjusted Sharpe ratio of the return distribution}
\usage{
AdjustedSharpeRatio(R, Rf = 0, FUN = "StdDev")
}
\arguments{
\item{R}{an xts, vector, matrix, data frame, timeSeries or zoo object of
asset returns}

\item{Rf}{the risk free rate}

\item{FUN}{one of "StdDev" or "VaR" or "ES" to use as the denominator for unadjusted Sharpe ratio, default="StdDev"}

}

\description{
Adjusted Sharpe ratio was introduced by Pezier and White (2006) to adjusts
for skewness and kurtosis by incorporating a penalty factor for negative skewness
and excess kurtosis.
}
\details{
\deqn{Adjusted Sharpe Ratio = SR * [1 + (\frac{S}{6}) * SR - (\frac{K - 3}{24}) * SR^2]}{Adjusted Sharpe ratio = SR x [1 + (S/6) x SR - ((K-3) / 24) x SR^2]}

where \eqn{SR} is the sharpe ratio with data annualized, \eqn{S} is the skewness and \eqn{K} is the kurtosis
}
\examples{
    data(assetReturns)

AdjustedSharpeRatio(assetReturns)

}

\references{
Carl Bacon, \emph{Practical portfolio performance measurement
and attribution}, second edition 2008 p.99.\cr

Pezier, Jaques and White, Anthony. 2006. The Relative Merits of Investable
Hedge Fund Indices and of Funds of Hedge Funds in Optimal Passive Portfolios.
Check \url{https://econpapers.repec.org/paper/rdgicmadp/icma-dp2006-10.htm} \cr
See also package \code{PerformanceAnalytics}.
}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}
