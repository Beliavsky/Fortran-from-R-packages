\name{varg}
\alias{varg}
\title{Value At Risk}
\description{Computes value at risk for a given ditribution}
\usage{
varg(spec, alpha, ...)
}
\arguments{
  \item{spec}{a character string specifying the distribution (for example, "norm" corresponds to the standard normal)}
  \item{alpha}{the probabilities associated with values at risk}
  \item{...}{other parameters}}
\value{An object of the same length as \code{alpha}, giving values at risk computed.}
\references{S. Chan and S. Nadarajah, Risk: An R package for risk measures, submitted}
\author{Stephen Chan, Saralees Nadarajah}
\examples{varg("norm", 0.9)}
\keyword{Value at risk}
