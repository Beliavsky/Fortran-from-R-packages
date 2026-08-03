\name{MartinRatio}
\alias{MartinRatio}

\encoding{latin1}

\title{Martin ratio of the return distribution}
\usage{
MartinRatio(R, Rf = 0)
}
\arguments{
\item{R}{an xts, vector, matrix, data frame, timeSeries or zoo object of
asset returns}

\item{Rf}{risk free rate, in same period as your returns}

}
\description{
To calculate Martin ratio we divide the difference of the portfolio return
and the risk free rate by the Ulcer index
}
\details{
\deqn{Martin ratio = \frac{r_P - r_F}{\sqrt{\sum^{n}_{i=1} \frac{{D'_i}^2}{n}}}}{Martin ratio = (rp - rf) / Ulcer index}

where \eqn{r_P} is the annualized portfolio return, \eqn{r_F} is the risk free
rate, \eqn{n} is the number of observations of the entire series, \eqn{D'_i} is
the drawdown since previous peak in period i
}
\examples{
  data(assetReturns)
	R=assetReturns[, -29]

# Not run
# MartinRatio(R)

}
\references{
Carl Bacon, \emph{Practical portfolio performance measurement
and attribution}, second edition 2008 p.91\cr
See also package \code{PerformanceAnalytics}.
}
\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}
