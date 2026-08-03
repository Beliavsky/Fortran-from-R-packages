\name{ProspectRatio}
\alias{ProspectRatio}

\encoding{latin1}
\title{Prospect ratio of the return distribution}
\usage{
ProspectRatio(R, MAR)
}
\arguments{
\item{R}{an xts, vector, matrix, data frame, timeSeries or zoo object of
asset returns}

\item{MAR}{the minimum acceptable return}

}
\description{
Prospect ratio is a ratio used to penalise loss since most people feel loss
greater than gain
}
\details{
\deqn{ProspectRatio(R) = \frac{\frac{1}{n}*\sum^{n}_{i=1}(Max(r_i,0)+2.25*Min(r_i,0) - MAR)}{\sigma_D}}{ProspectRatio(R) = (1/n * sum(Max(ri,0) + 2.25 * Min(ri,0)) - MAR) / DownsideRisk}

where \eqn{n} is the number of observations of the entire series, MAR is the minimum acceptable return and \eqn{\sigma_D} is the downside risk
}
\examples{
  data(assetReturns)
	R=assetReturns[, -29]

  ProspectRatio(R, MAR=0)

}
\references{
Carl Bacon, \emph{Practical portfolio performance measurement
and attribution}, second edition 2008 p.100\cr
See also package \code{PerformanceAnalytics}.
}
\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}
