\name{luceg4}
\alias{luceg4}
\title{Luce (1980)'s Fourth Risk Measure}
\description{Computes the fourth risk measure due to Luce (1980)}
\usage{
luceg4(spec, a, b, aa, bb, ...)
}
\arguments{
  \item{spec}{a character string specifying the distribution (for example, "norm" corresponds to the standard normal)}
  \item{aa}{a positive valued parameter, see Chan and Nadarajah for details}
  \item{bb}{a positive valued parameter, see Chan and Nadarajah for details}
  \item{a}{the lower end point of the distribution specified by \code{spec}}
  \item{b}{the upper end point of the distribution specified by \code{spec}}
  \item{...}{other parameters}}
\value{A scalar, giving Luce (1980)'s fourth risk measure of the distribution specified by \code{spec}}
\references{S. Chan and S. Nadarajah, Risk: An R package for risk measures, submitted

R. D. Luce, Several possible measures of risk, Theory and Decision, 12, 1980, 217-228
<DOI:10.1007/BF00135033>}
\author{Stephen Chan, Saralees Nadarajah}
\examples{luceg4("norm",-Inf, Inf, 1, 0)}
\keyword{Luce (1980)'s fourth risk measure}
