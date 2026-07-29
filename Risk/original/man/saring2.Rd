\name{saring2}
\alias{saring2}
\title{Sarin (1987)'s Second Risk Measure}
\description{Computes the second risk measure due to Sarin (1987)}
\usage{
saring2(spec, a, b, aa, bb1, bb2, ...)
}
\arguments{
  \item{spec}{a character string specifying the distribution (for example, "norm" corresponds to the standard normal)}
  \item{aa}{a positive real valued parameter, see Chan and Nadarajah for details}
  \item{bb1}{a positive real valued parameter, see Chan and Nadarajah for details}
  \item{bb2}{a positive real valued parameter, see Chan and Nadarajah for details}
  \item{a}{the lower end point of the distribution specified by \code{spec}}
  \item{b}{the upper end point of the distribution specified by \code{spec}}
  \item{...}{other parameters}}
\value{A scalar, giving Sarin (1987)'s second risk measure of the distribution specified by \code{spec}}
\references{S. Chan and S. Nadarajah, Risk: An R package for risk measures, submitted

R. K. Sarin, Some extensions of Luce's measures of risk, Theory and Decision, 22, 1987, 125-141
<DOI:10.1007/BF00126387>}
\author{Stephen Chan, Saralees Nadarajah}
\examples{saring2("norm",-Inf, Inf, 1, 1, 1)}
\keyword{Sarin (1987)'s second risk measure}
