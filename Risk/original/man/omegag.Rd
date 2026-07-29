\name{omegag}
\alias{omegag}
\title{Omega Risk Measure Due To Shadwick And Keating (2002)}
\description{Computes the omega risk measure for a given ditribution}
\usage{
omegag(spec, alpha, a, b, ...)
}
\arguments{
  \item{spec}{a character string specifying the distribution (for example, "norm" corresponds to the standard normal)}
  \item{alpha}{a real valued parameter, see Chan and Nadarajah for details}
  \item{a}{the lower end point of the distribution specified by \code{spec}}
  \item{b}{the upper end point of the distribution specified by \code{spec}}
  \item{...}{other parameters}}
\value{An object of the same length as \code{alpha}, giving the omega risk measure of the distribution specified by \code{spec}}
\references{S. Chan and S. Nadarajah, Risk: An R package for risk measures, submitted

W. F. Shadwick and C. Keating,
A universal performance measure,
Journal of  Performance Measurement, 2002}
\author{Stephen Chan, Saralees Nadarajah}
\examples{omegag("norm", 2, -Inf, Inf)}
\keyword{Omega risk measure}
