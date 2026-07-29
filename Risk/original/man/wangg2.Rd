\name{wangg2}
\alias{wangg2}
\title{Wang (1998)'s Second Risk Measure}
\description{Computes the second risk measure due to Wang (1998)}
\usage{
wangg2(spec, alpha, a, b, ...)
}
\arguments{
  \item{spec}{a character string specifying the distribution (for example, "norm" corresponds to the standard normal)}
  \item{alpha}{a real valued parameter taking values in (0, 1), see Chan and Nadarajah for details}
  \item{a}{the lower end point of the distribution specified by \code{spec}}
  \item{b}{the upper end point of the distribution specified by \code{spec}}
  \item{...}{other parameters}}
\value{An object of the same length as \code{alpha}, giving Wang (1998)'s second risk measure of the distribution specified by \code{spec}}
\references{S. Chan and S. Nadarajah, Risk: An R package for risk measures, submitted

S. Wang, An actuarial index of the right-tail risk,
North American Actuarial Journal, 2, 1998, 88-101
<DOI:10.1080/10920277.1998.10595708>}
\author{Stephen Chan, Saralees Nadarajah}
\examples{wangg2("lnorm", 0.9, 0, Inf)}
\keyword{Wang (1998)'s second risk measure}
