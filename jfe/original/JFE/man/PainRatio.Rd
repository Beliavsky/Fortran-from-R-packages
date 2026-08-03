\name{PainRatio}
\alias{PainRatio}
\encoding{latin1}

\title{Pain ratio of the return distribution}
\usage{
PainRatio(R, Rf = 0)
}
\arguments{
\item{R}{an xts, vector, matrix, data frame, timeSeries or zoo object of
asset returns}

\item{Rf}{risk free rate, in same period as your returns}

}
\description{
To calculate Pain ratio we divide the difference of the portfolio return
and the risk free rate by the Pain index
}
\details{
\deqn{Pain ratio = \frac{r_P - r_F}{\sum^{n}_{i=1} \frac{\mid D'_i \mid}{n}}}{Pain ratio = (rp - rf) / Pain index}

where \eqn{r_P} is the annualized portfolio return, \eqn{r_F} is the risk free
rate, \eqn{n} is the number of observations of the entire series, \eqn{D'_i} is
the drawdown since previous peak in period i
}
\examples{
  data(assetReturns)
	assetReturns=assetReturns["2011::2018"] #short sample for fast example
	R=assetReturns[, -29]
  PainRatio(R)

}
\references{
Carl Bacon, \emph{Practical portfolio performance measurement
and attribution}, second edition 2008 p.91\cr
See also package \code{PerformanceAnalytics}.
}

\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}
