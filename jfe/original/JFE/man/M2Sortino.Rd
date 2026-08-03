\name{M2Sortino}
\alias{M2Sortino}

\encoding{latin1}

\title{M squared for Sortino of the return distribution}
\usage{
M2Sortino(Ra, Rb, MAR = 0)
}
\arguments{
\item{Ra}{an xts, vector, matrix, data frame, timeSeries or zoo object of
asset return}

\item{Rb}{return vector of the benchmark asset}

\item{MAR}{the minimum acceptable return}

}
\description{
M squared for Sortino is a M^2 calculated for Downside risk instead of Total Risk
}
\details{
\deqn{M^2_S = r_P + Sortino ratio * (\sigma_{DM} - \sigma_D)}{M^2 (Sortino) = Rp + Sortino ratio * (DownsideRiskBenchmark - DownsideRiskPortfolio)}

where \eqn{M^2_S} is MSquared for Sortino, \eqn{r_P} is the annualised portfolio return,
\eqn{\sigma_{DM}} is the benchmark annualised downside risk and \eqn{D} is the portfolio
annualised downside risk
}
\examples{

  data(assetReturns)
	Ra=assetReturns[, -29]
	Rb=assetReturns[,29] #DJI

M2Sortino(Ra, Rb, MAR=0)

}
\references{
Carl Bacon, \emph{Practical portfolio performance measurement
and attribution}, second edition 2008 p.102-103\cr
See aslo package \code{PerformanceAnalytics}.
}
\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

