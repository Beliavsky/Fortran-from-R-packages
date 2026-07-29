\name{tcm}
\alias{tcm}
\title{Tail Conditional Mean Due To Kou et al. (2013)}
\description{Computes tail conditional median for a given ditribution}
\usage{
tcm(spec, alpha, ...)
}
\arguments{
  \item{spec}{a character string specifying the distribution (for example, "norm" corresponds to the standard normal)}
  \item{alpha}{the probabilities associated with tail conditional median}
  \item{...}{other parameters}}
\value{An object of the same length as \code{alpha}, giving tail conditional medians computed.}
\references{S. Chan and S. Nadarajah, Risk: An R package for risk measures, submitted

S. Kou, X. Peng and C. C. Heyde,
External risk measures and Basel accords, Mathematics of Operations Research, 38, 2013, 393-417
<DOI:10.1287/moor.1120.0577>}
\author{Stephen Chan, Saralees Nadarajah}
\examples{tcm("norm", 0.9)}
\keyword{Tail conditional median}
