\name{rebuild.cov}
\alias{rebuild.cov}
\alias{rebuild.invcov}
\alias{decompose.cov}
\alias{decompose.invcov}


\title{Rebuild and Decompose the (Inverse) Covariance Matrix}
\usage{

rebuild.cov(r, v)
rebuild.invcov(pr, pv)
decompose.cov(m)
decompose.invcov(m)
}
\arguments{
  \item{r}{correlation matrix}
  \item{v}{variance vector}
  \item{pr}{partial correlation matrix}
  \item{pv}{partial variance vector}
  \item{m}{a covariance or a concentration matrix}
  
}
\description{
  
  \code{rebuild.cov} takes a correlation matrix and a vector with variances
  and reconstructs the corresponding covariance matrix.
  
  Conversely, \code{decompose.cov} decomposes 
  a covariance matrix into correlations and variances.
  
  \code{decompose.invcov} decomposes a concentration matrix (=inverse covariance
  matrix) into partial correlations and partial variances.  
  
  \code{rebuild.invcov} takes a partial correlation matrix and a vector 
  with partial variances and reconstructs the corresponding concentration matrix.
}

\details{
  The diagonal elements of the concentration matrix (=inverse covariance matrix)
  are the precisions, and the off-diagonal elements
  are the concentrations. Thus, the partial variances
  correspond to the inverse precisions, and the partial correlations to the
  negative standardized concentrations.
}


\value{

  \code{rebuild.cov} and \code{rebuild.invcov} return a matrix.
  
   \code{decompose.cov} and \code{decompose.invcov} return a list containing
   a matrix and a vector.
}


\author{
  Korbinian Strimmer (\url{https://strimmerlab.github.io}).
}
\seealso{\code{\link{cor}}, \code{\link{cov}}, \code{\link{pcor.shrink}}}

\examples{
# load corpcor library
library("corpcor")

# a correlation matrix and some variances
r = matrix(c(1, 1/2, 1/2, 1),  nrow = 2, ncol=2)
r
v = c(2, 3)

# construct the associated covariance matrix
c = rebuild.cov(r, v)
c

# decompose into correlations and variances
decompose.cov(c)


# the corresponding concentration matrix
conc = pseudoinverse(c) 
conc

# decompose into partial correlation matrix and partial variances
tmp = decompose.invcov(conc)
tmp
# note: because this is an example with two variables,
# the partial and standard correlations are identical!


# reconstruct the concentration matrix from partial correlations and
# partial variances 
rebuild.invcov(tmp$pr, tmp$pv)

}
\keyword{multivariate}
