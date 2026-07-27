\name{riskfun}
\alias{riskfun}
\title{
Portfolio Risk Measures
}
\description{
Calculates a given portfolio risk/deviation measure given a set of weights
and matrix of returns, possible representing a forecast scenario.
}
\usage{
riskfun(weights, Data, risk = c("mad", "ev", "minimax", "cvar", "cdar", "lpm"), 
benchmark = NULL, alpha = 0.05, moment = 1, threshold = 0, VaR = NULL, DaR = NULL)
}
\arguments{
  \item{weights}{ vector of weights.}
  \item{Data}{ Matrix of returns.}
  \item{risk}{ Choice of measure.}
  \item{benchmark}{ (Optional) vector of benchmark returns with same number of
  rows as Data.}
  \item{alpha}{ The lower quantile for the \dQuote{cvar} and \dQuote{cdar} 
  measures.}
  \item{moment}{ The \dQuote{lpm} measure moment.}
  \item{threshold}{ The \dQuote{lpm} measure threshold. A value of 999 will 
  subtract the portfolio mean.}
  \item{VaR}{ (Optional) The pre-calculated VaR for the \dQuote{cvar} measure.}
  \item{DaR}{ (Optional) The pre-calculated DaR for the \dQuote{cdar} measure.}
}
\details{
A simple utility function for the calculation and understanding of some of the 
risk and deviation measures implemented in the package.
}
\value{
A numeric value representing the risk/deviation measure.
}
\author{
Alexios Galanos
}
\keyword{optimize}
