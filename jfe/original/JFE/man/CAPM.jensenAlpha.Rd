\name{CAPM.jensenAlpha}
\alias{CAPM.jensenAlpha}

\encoding{latin1}

\title{Jensen's alpha of the return distribution}
\usage{
CAPM.jensenAlpha(Ra, Rb, Rf = 0)
}
\arguments{
\item{Ra}{an xts, vector, matrix, data frame, timeSeries or zoo object of
asset returns}

\item{Rb}{return vector of the benchmark asset}

\item{Rf}{risk free rate, in same period as your returns}

}
\description{
The Jensen's alpha is the intercept of the regression equation in the Capital
Asset Pricing Model and is in effect the exess return adjusted for systematic risk.
}
\details{
\deqn{\alpha = r_p - r_f - \beta_p * (b - r_f)}{alpha = r_p - r_f - beta_p * (b - r_f)}

where \eqn{r_f} is the risk free rate, \eqn{\beta_r} is the regression beta,
\eqn{r_p} is the portfolio return and b is the benchmark return
}
\examples{

  data(assetReturns)
	assetReturns=assetReturns["2011::2018"] #short sample for fast example
	Ra=assetReturns[, -29]
	Rb=assetReturns[,29] #DJI

CAPM.jensenAlpha(Ra, Rb)

}
\references{
Carl Bacon, \emph{Practical portfolio performance measurement
and attribution}, second edition 2008 p.72\cr
See also package \code{PerformanceAnalytics}.
}
\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}
