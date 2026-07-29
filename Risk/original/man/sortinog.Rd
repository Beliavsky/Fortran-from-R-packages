\name{sortinog}
\alias{sortinog}
\title{Sortino Ratio Due To Rollinger And Hoffman (2013)}
\description{Computes the Sortino ratio for a given ditribution}
\usage{
sortinog(spec, alpha, a, b, ...)
}
\arguments{
  \item{spec}{a character string specifying the distribution (for example, "norm" corresponds to the standard normal)}
  \item{alpha}{a real valued parameter, see Chan and Nadarajah for details}
  \item{a}{the lower end point of the distribution specified by \code{spec}}
  \item{b}{the upper end point of the distribution specified by \code{spec}}
  \item{...}{other parameters}}
\value{An object of the same length as \code{alpha}, giving the Sortino ratio of the distribution specified by \code{spec}}
\references{S. Chan and S. Nadarajah, Risk: An R package for risk measures, submitted

T. Rollinger and S. Hoffman,
Sortino ratio: A better measure of risk, Risk Management, 40-42, 2013}
\author{Stephen Chan, Saralees Nadarajah}
\examples{sortinog("norm", 2, -Inf, Inf)}
\keyword{Sortino ratio}
