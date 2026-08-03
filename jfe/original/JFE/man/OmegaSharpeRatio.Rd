\name{OmegaSharpeRatio}
\alias{OmegaSharpeRatio}

\encoding{latin1}

\title{Omega-Sharpe ratio of the return distribution}
\usage{
OmegaSharpeRatio(R, MAR = 0)
}
\arguments{
\item{R}{an xts, vector, matrix, data frame, timeSeries or zoo object of
asset returns}

\item{MAR}{Minimum Acceptable Return, in the same periodicity as your
returns}

}

\description{
The Omega-Sharpe ratio is a conversion of the omega ratio to a ranking statistic
in familiar form to the Sharpe ratio.
}
\details{
To calculate the Omega-Sharpe ration we subtract the target (or Minimum
Acceptable Returns (MAR)) return from the portfolio return and we divide
it by the opposite of the Downside Deviation.

\deqn{OmegaSharpeRatio(R,MAR) = \frac{r_p - r_t}{\sum^n_{t=1}\frac{max(r_t - r_i, 0)}{n}}}{OmegaSharpeRatio(R,MAR) = (Rp - Rt) / -DownsidePotential(R,MAR)}

where \eqn{n} is the number of observations of the entire series
}
\examples{

  data(assetReturns)
	R=assetReturns[, -29]
  OmegaSharpeRatio(R)

}
\references{
Carl Bacon, \emph{Practical portfolio performance measurement
and attribution}, second edition 2008, p.95\cr
See also package \code{PerformanceAnalytics}.
}
\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

