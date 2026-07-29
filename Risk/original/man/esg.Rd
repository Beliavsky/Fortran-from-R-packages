\name{esg}
\alias{esg}
\title{Expected Shortfall Due To Artzner et al. (1999)}
\description{Computes expected shortfall for a given ditribution}
\usage{
esg(spec, alpha, ...)
}
\arguments{
  \item{spec}{a character string specifying the distribution (for example, "norm" corresponds to the standard normal)}
  \item{alpha}{the probabilities associated with expected shortfall}
  \item{...}{other parameters}}
\value{An object of the same length as \code{alpha}, giving expected shortfall computed.}
\references{S. Chan and S. Nadarajah, Risk: An R package for risk measures, submitted

P. Artzner, F. Delbaen, J. M. Eber and D. Heath,
Coherent measures of risk, Mathematical Finance, 9, 1999, 203-228
<DOI:10.1111/1467-9965.00068>}
\author{Stephen Chan, Saralees Nadarajah}
\examples{esg("norm", 0.9)}
\keyword{Expected shortfall}
