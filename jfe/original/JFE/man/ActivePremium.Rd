\name{ActivePremium}
\alias{ActivePremium}

\encoding{latin1}

\title{Active Premium or Active Return}

\description{
The return on an investment's annualized return minus the benchmark's
annualized return.
}

\usage{
ActivePremium(Ra, Rb, scale = NA)
}
\arguments{
\item{Ra}{return vector of the portfolio}

\item{Rb}{return vector of the benchmark asset}

\item{scale}{number of periods in a year
(daily scale = 252, monthly scale = 12, quarterly scale = 4)}

}

\details{
Active Premium = Investment's annualized return - Benchmark's annualized
return.
With a view to speeding computation. I re-write the code of some ratios of the package PerformanceAnalytics,  and use the same name for comparing the performance enhancing. Interested readers may compare speed improvement with the use of system.time().
}
\examples{

  data(assetReturns)
	assetReturns=assetReturns["2015::2018"] #short sample for fast example
#	Ra=assetReturns[, -29]
#	Rb=assetReturns[,29] #DJI
#  ActivePremium(Ra, Rb)

}
\references{
Sharpe, W.F. The Sharpe Ratio,\emph{Journal of Portfolio Management}, Fall 1994, 49-58.\cr
See aslo package \code{PerformanceAnalytics}.
}
\seealso{
\code{\link{Return.annualized}}, \cr
\code{\link{InformationRatio}}, \cr
\code{\link{TrackingError}}, \cr
}


\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}
