\name{table.AnnualizedReturns}
\alias{table.AnnualizedReturns}
\encoding{latin1}
\title{Annualized Returns Summary: Statistics and Stylized Facts}
\usage{
table.AnnualizedReturns(R, scale = NA, Rf = 0, geometric = TRUE,
  digits = 4)
}
\arguments{
\item{R}{an xts, vector, matrix, data frame, timeSeries or zoo object of
asset returns}

\item{scale}{number of periods in a year (daily scale = 252, monthly scale =
12, quarterly scale = 4)}

\item{Rf}{risk free rate, in same period as your returns}

\item{geometric}{utilize geometric chaining (TRUE) or simple/arithmetic chaining (FALSE) to aggregate returns,
default TRUE}

\item{digits}{number of digits to round results to}
}
\description{
Table of Annualized Return, Annualized Std Dev, and Annualized Sharpe
}
\examples{

  data(assetReturns)
	Ra=assetReturns[, -29]
  table.AnnualizedReturns(R=Ra)

}
\seealso{
\code{\link{Return.annualized}} \cr
\code{\link{SharpeRatio.annualized}}
}
\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

