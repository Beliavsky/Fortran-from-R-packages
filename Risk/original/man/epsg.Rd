\name{epsg}
\alias{epsg}
\title{Expected Proportional Shortfall Due To Belzunce et al. (2012)}
\description{Computes expected proportional shortfall for a given ditribution}
\usage{
epsg(spec, alpha, ...)
}
\arguments{
  \item{spec}{a character string specifying the distribution (for example, "norm" corresponds to the standard normal)}
  \item{alpha}{the probabilities associated with expected proportional shortfalls}
  \item{...}{other parameters}}
\value{An object of the same length as \code{alpha}, giving expected proportional shortfalls computed.}
\references{S. Chan and S. Nadarajah, Risk: An R package for risk measures, submitted

F. Belzunce, J. F. Pinar, J. M. Ruiz and M. A. Sordo,
Comparison of risks based on the expected proportional shortfall,
Insurance: Mathematics and Economics, 51, 2012,  292-302
<DOI:10.1016/j.insmatheco.2012.05.003>}
\author{Stephen Chan, Saralees Nadarajah}
\examples{epsg("norm", 0.9)}
\keyword{Expected proportional shortfall}
