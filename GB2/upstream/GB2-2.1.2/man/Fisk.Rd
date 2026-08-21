\name{Fisk}
\Rdversion{2.1}
\alias{Fisk}
\alias{fisk}
\alias{fiskh}
\title{
Parameters of the Fisk Distribution
}
\description{
Calculation of the parameters \eqn{a} and \eqn{b} of the Fisk distribution, which is a GB2 distribution with \eqn{p = q = 1}. If \eqn{m} and \eqn{v} denote, respectively, the mean and variance of \eqn{log(z)}, 
then \eqn{\hat{a} = \pi/\sqrt{3*v}} and \eqn{\hat{b} = \exp(m)}.  
}
\usage{
fisk(z, w=rep(1, length(z)))
fiskh(z, w=rep(1, length(z)), hs=rep(1, length(z)))
}
\arguments{
  \item{z}{numeric; vector of data values.}
  \item{w}{numeric; vector of weights. Must have the same length as \code{z}. By default \code{w} is a vector of 1.}
  \item{hs}{numeric; vector of household sizes. Must have the same length as \code{z}. By default \code{hs} is a vector of 1.}
}
\details{
Function \code{fisk} first calculates the mean and variance of \eqn{log(z)} and next the values of \eqn{a} and \eqn{b} under the Fisk distribution.  
Function \code{fiskh} first calculates the mean and variance of \eqn{log(z)}, assuming a sample of households, and next the values of \eqn{a} and \eqn{b} under the Fisk distribution. 
}
\value{\code{fisk} and \code{fiskh} return vectors of length 4 containing the estimated parameters \eqn{a}, \eqn{b}, as well as \eqn{p=1} and \eqn{q=1}. 
}
\references{
 Graf, M., Nedyalkova, D., Muennich, R., Seger, J. and Zins, S. (2011)
 AMELI Deliverable 2.1: Parametric Estimation of Income Distributions and
 Indicators of Poverty and Social Exclusion.
 \emph{Technical report}, AMELI-Project.
}
\seealso{\code{\link{optim}} for the general-purpose optimization
}
\author{
Monique Graf
}
\examples{
library(laeken)
data(eusilc)

# Income
inc <- as.vector(eusilc$eqIncome)

# Weights
w <- eusilc$rb050

#Fisk parameters
fpar <- fisk(inc, w)
}
\keyword{distribution}

