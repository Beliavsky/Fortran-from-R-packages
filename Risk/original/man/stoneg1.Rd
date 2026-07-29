\name{stoneg1}
\alias{stoneg1}
\title{Stone (1973)'s First Risk Measure}
\description{Computes the first risk measure due to Stone (1973)}
\usage{
stoneg1(spec, x0, k, a, b, ...)
}
\arguments{
  \item{spec}{a character string specifying the distribution (for example, "norm" corresponds to the standard normal)}
  \item{x0}{a real valued parameter, see Chan and Nadarajah for details}
  \item{k}{a positive valued parameter, see Chan and Nadarajah for details}
  \item{a}{the lower end point of the distribution specified by \code{spec}}
  \item{b}{the upper end point of the distribution specified by \code{spec}}
  \item{...}{other parameters}}
\value{A scalar, giving Stone (1973)'s first risk measure of the distribution specified by \code{spec}}
\references{S. Chan and S. Nadarajah, Risk: An R package for risk measures, submitted

B. K. Stone,
A general class of three-parameter risk measuresm, The Journal of Finance, 28, 1973, 675-685
<DOI:10.2307/2978638>}
\author{Stephen Chan, Saralees Nadarajah}
\examples{stoneg1("norm", 8, 3, -Inf, Inf)}
\keyword{Stone (1973)'s first risk measure}
