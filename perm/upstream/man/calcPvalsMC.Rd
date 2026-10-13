\name{calcPvalsMC}
\alias{calcPvalsMC}
\alias{calcPvals}
\alias{pCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Internal functions: Calculate p-values from Complete Enumeration or Monte Carlo}
\description{
These internal functions are not meant to be called directly. This help is to explain the 
functions. The function \code{calcPvals} calculates p-values from a set of either 
completely enumerated test statistics (in which case all wgts are equal to one) or 
a set of test statistics from a network algorithm (in which case the wgts come from the network algorithm).
The function \code{calcPvalsMC} calculates p-values and the associated confidence intervals 
from a Monte Carlo simulation. The p-values are estimated by (S+1)/(N+1) where S is the number of 
more extreme test statistics, and N is the number of Monte Carlo simulations. This estimator is 
guaranteed to be valid,  and the confidence 
interval is the appropriate one associated with this estimator (see e.g., Fay, Kim, and Hachey, 2007). 
}
\usage{
calcPvalsMC(Tj, T0, digits, alternative = NULL, twosidedTstat = FALSE, p.conf.level=0.99)
calcPvals(Tj, T0, digits, denom = NULL, wgts = NULL, twosidedTstat = FALSE)
pCI(S, N, p.conf.level = 0.99)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Tj}{ numeric vector of values of test statistic }
  \item{T0}{ observed value of test statistic}
  \item{digits}{ number of significant digits for rounding (see details)}
  \item{alternative}{character vector: one of "less", "greater", "two.sided"}
  \item{twosidedTstat}{logical, TRUE denotes test statistic is twosided}
  \item{denom}{total number of test statistics}
  \item{wgts}{vector of weights the same length as Tj, NULL gives all ones}
  \item{S}{number of Monte Carlo replications equal or more extreme than T0}
\item{N}{number of Monte Carlo replications}
\item{p.conf.level}{confidence level for confidence intervals about Monte Carlo p-values}
}
%\details{  ~~ If necessary, more details than the description above ~~}
%\value{ }
\references{Fay, MP, Kim, H-J, Hachey, M (2007). 
Journal of Computational and Graphical Statistics.
16:946-957. Available at
\url{https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2467508/}
}
%\author{ ~~who you are~~ }
%\note{ ~~further notes~~ }
\seealso{Called from \code{\link{permTS}},  \code{\link{permKS}},  \code{\link{permTREND}}  }
%\examples{
\keyword{internal}
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
