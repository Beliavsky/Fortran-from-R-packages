\name{TreynorRatio}
\alias{TreynorRatio}
\encoding{latin1}
\title{calculate Treynor Ratio or modified Treynor Ratio of excess return over CAPM beta}
\usage{
TreynorRatio(Ra, Rb, Rf = 0, scale = NA, modified = FALSE)
}
\arguments{
\item{Ra}{an xts, vector, matrix, data frame, timeSeries or zoo object of
asset returns}

\item{Rb}{return vector of the benchmark asset}

\item{Rf}{risk free rate, in same period as your returns}

\item{scale}{number of periods in a year (daily scale = 252, monthly scale =
12, quarterly scale = 4)}

\item{modified}{a boolean to decide whether to return the Treynor ratio or
Modified Treynor ratio}
}
\description{
The Treynor ratio is similar to the Sharpe Ratio, except it uses beta as the
volatility measure (to divide the investment's excess return over the beta).
}
\details{
To calculate modified Treynor ratio, we divide the numerator by the systematic risk
instead of the beta.

Equation:
\deqn{TreynorRatio = \frac{\overline{(R_{a}-R_{f})}}{\beta_{a,b}}}{(mean(Ra-Rf))/(Beta(Ra,Rb))}
\deqn{ModifiedTreynorRatio = \frac{r_p - r_f}{\sigma_s}}{ModifiedTreynorRatio = (Rp - Rf)/sytematic risk}
}
\examples{

  data(assetReturns)
	assetReturns=assetReturns["2011::2018"] #short sample for fast example
	Ra=assetReturns[, -29]
	Rb=assetReturns[,29] #DJI

  TreynorRatio(Ra, Rb)

}
\references{
\url{https://en.wikipedia.org/wiki/Treynor_ratio},
Carl Bacon, \emph{Practical portfolio performance measurement
and attribution}, second edition 2008 p.77\cr
See also package \code{PerformanceAnalytics}.
}
\seealso{
\code{\link{SharpeRatio}} \code{\link{SortinoRatio}}

}
\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

