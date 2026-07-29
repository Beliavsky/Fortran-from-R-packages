\name{bvar}
\alias{bvar}
\title{Beyond Value At Risk Due To Longin (2001)}
\description{Computes beyond value at risk for a given ditribution}
\usage{
bvar(spec, alpha, a, ...)
}
\arguments{
  \item{spec}{a character string specifying the distribution (for example, "norm" corresponds to the standard normal)}
  \item{alpha}{the probabilities associated with beyon values at risk}
  \item{a}{the lower end point of the distribution specified by \code{spec}}
  \item{...}{other parameters}}
\value{An object of the same length as \code{alpha}, giving beyond values ar risk computed.}
\references{S. Chan and S. Nadarajah, Risk: An R package for risk measures, submitted

F. M. Longin,
Beyond the VaR,
Journal of Derivatives, 8, 2001, 36-48
<DOI:10.3905/jod.2001.319161>}
\author{Stephen Chan, Saralees Nadarajah}
\examples{bvar("norm", 0.9, a=-Inf)}
\keyword{Beyond value ar risk}
