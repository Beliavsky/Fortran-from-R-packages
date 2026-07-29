\name{expect}
\alias{expect}
\title{Expectation}
\description{Computes expectation for a given ditribution}
\usage{
expect(spec, a, b, ...)
}
\arguments{
  \item{spec}{a character string specifying the distribution (for example, "norm" corresponds to the standard normal)}
  \item{a}{the lower end point of the distribution specified by \code{spec}}
  \item{b}{the upper end point of the distribution specified by \code{spec}}
  \item{...}{other parameters}}
\value{A scalar, giving the expected value of the distribution specified by \code{spec}}
\references{S. Chan and S. Nadarajah, Risk: An R package for risk measures, submitted}
\author{Stephen Chan, Saralees Nadarajah}
\examples{expect("norm", -Inf, Inf)}
\keyword{Expectation}
