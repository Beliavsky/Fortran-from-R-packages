\name{BernardoLedoitRatio}
\alias{BernardoLedoitRatio}

\encoding{latin1}

\title{Bernardo and Ledoit ratio of the return distribution}
\usage{
BernardoLedoitRatio(R)
}
\arguments{
\item{R}{an xts, vector, matrix, data frame, timeSeries or zoo object of
asset returns}

}
\description{
To calculate Bernardo and Ledoit ratio we take the sum of the subset of
returns that are above 0 and we divide it by the opposite of the sum of
the subset of returns that are below 0
}
\details{
\deqn{BernardoLedoitRatio(R) = \frac{\frac{1}{n}\sum^{n}_{t=1}{max(R_{t},0)}}{\frac{1}{n}\sum^{n}_{t=1}{max(-R_{t},0)}}}{BernardoLedoitRatio(R) = 1/n*sum(t=1..n)(max(R(t),0)) / 1/n*sum(t=1..n)(max(-R(t),0))}

where \eqn{n} is the number of observations of the entire series
}
\examples{
    data(assetReturns)
BernardoLedoitRatio(R=assetReturns)



}
\references{
Carl Bacon, \emph{Practical portfolio performance measurement
and attribution}, second edition 2008 p.95.\cr
See also package \code{PerformanceAnalytics}.
}
\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}
