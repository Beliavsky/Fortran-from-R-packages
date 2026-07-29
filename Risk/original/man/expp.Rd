\name{expp}
\alias{expp}
\title{Expectiles Due To Newey And Powell (1987)}
\description{Computes expectiles for a given ditribution}
\usage{
expp(spec, alpha, a, b, ...)
}
\arguments{
  \item{spec}{a character string specifying the distribution (for example, "norm" corresponds to the standard normal)}
  \item{alpha}{the probabilities associated with expectiles}
  \item{a}{the lower end point of the distribution specified by \code{spec}}
  \item{b}{the upper end point of the distribution specified by \code{spec}}
  \item{...}{other parameters}}
\value{An object of the same length as \code{alpha}, giving expectiles computed.}
\references{S. Chan and S. Nadarajah, Risk: An R package for risk measures, submitted

W. K. Newey and J. L. Powell,
Asymmetric least squares estimation and testing. Econometrica, 55, 1987, 819-847
<DOI:10.2307/1911031>}
\author{Stephen Chan, Saralees Nadarajah}
\examples{expp("norm", 0.9, a=-Inf, b=Inf)}
\keyword{Expectiles}
