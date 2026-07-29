\name{expvar}
\alias{expvar}
\title{An Elementary Risk Measure Due To Ahmadi-Javid (2012)}
\description{Computes the elementary risk measure for a given ditribution}
\usage{
expvar(spec, alpha, a, b, ...)
}
\arguments{
  \item{spec}{a character string specifying the distribution (for example, "norm" corresponds to the standard normal)}
  \item{alpha}{a positive valued parameter, see Chan and Nadarajah for details}
  \item{a}{the lower end point of the distribution specified by \code{spec}}
  \item{b}{the upper end point of the distribution specified by \code{spec}}
  \item{...}{other parameters}}
\value{An object of the same length as \code{alpha}, giving the elementary risk measure of the distribution specified by \code{spec}}
\references{S. Chan and S. Nadarajah, Risk: An R package for risk measures, submitted

A. Ahmadi-Javid,
Entropic value-at-risk: A new coherent risk measure.
Journal of Optimization Theory and Applications, 155, 2012, 1105-1123
<DOI:10.1007/s10957-011-9968-2>}
\author{Stephen Chan, Saralees Nadarajah}
\examples{expvar("norm", 0.9, -Inf, Inf)}
\keyword{Elementary risk measure}
