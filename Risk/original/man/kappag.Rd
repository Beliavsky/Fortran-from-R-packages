\name{kappag}
\alias{kappag}
\title{Kappa Risk Measure Due To Kaplan And Knowles (2004)}
\description{Computes the Kappa risk measure for a given ditribution}
\usage{
kappag(spec, alpha, n, a, b, ...)
}
\arguments{
  \item{spec}{a character string specifying the distribution (for example, "norm" corresponds to the standard normal)}
  \item{alpha}{a real valued parameter, see Chan and Nadarajah for details}
  \item{n}{a positive integer valued parameter, see Chan and Nadarajah for details}
  \item{a}{the lower end point of the distribution specified by \code{spec}}
  \item{b}{the upper end point of the distribution specified by \code{spec}}
  \item{...}{other parameters}}
\value{An object of the same length as \code{alpha}, giving the Kappa risk measure of the distribution specified by \code{spec}}
\references{S. Chan and S. Nadarajah, Risk: An R package for risk measures, submitted

P. D. Kaplan and J. A. Knowles,
Kappa: A generalized downside risk-adjusted performance measure,
Miscellaneous Publication, Morningstar Associates and York Hedge Fund Strategies, 2004}
\author{Stephen Chan, Saralees Nadarajah}
\examples{kappag("norm", 2, 5, -Inf, Inf)}
\keyword{Kappa risk measure}
