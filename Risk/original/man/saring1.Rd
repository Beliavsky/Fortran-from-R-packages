\name{saring1}
\alias{saring1}
\title{Sarin (1987)'s First Risk Measure}
\description{Computes the first risk measure due to Sarin (1987)}
\usage{
saring1(spec, a, b, k, c, ...)
}
\arguments{
  \item{spec}{a character string specifying the distribution (for example, "norm" corresponds to the standard normal)}
  \item{k}{a non-zero real valued parameter, see Chan and Nadarajah for details}
  \item{c}{a non-zero real valued parameter, see Chan and Nadarajah for details}
  \item{a}{the lower end point of the distribution specified by \code{spec}}
  \item{b}{the upper end point of the distribution specified by \code{spec}}
  \item{...}{other parameters}}
\value{A scalar, giving Sarin (1987)'s first risk measure of the distribution specified by \code{spec}}
\references{S. Chan and S. Nadarajah, Risk: An R package for risk measures, submitted

R. K. Sarin, Some extensions of Luce's measures of risk, Theory and Decision, 22, 1987, 125-141
<DOI:10.1007/BF00126387>}
\author{Stephen Chan, Saralees Nadarajah}
\examples{saring1("norm", -Inf, Inf, 1, 0)}
\keyword{Sarin (1987)'s first risk measure}
