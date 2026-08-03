\name{wt.scale}
\alias{wt.scale}
\alias{wt.var}
\alias{wt.moments}

\title{Weighted Expectations and Variances}

\description{
  \code{wt.var} estimate the unbiased variance taking into account 
  data weights.

  \code{wt.moments} produces the weighted mean and weighted variance
  for each column of a matrix.
  
  \code{wt.scale} centers and standardized a matrix using 
  the weighted means and variances.
 }

\usage{
  wt.var(xvec, w) 
  wt.moments(x, w)
  wt.scale(x, w, center=TRUE, scale=TRUE)
}

\arguments{
  \item{xvec}{a vector}
  \item{x}{a matrix}
  \item{w}{data weights}
  \item{center}{logical value}
  \item{scale}{logical value}
}


\value{
  A rescaled matrix (\code{wt.scale}), a list containing the column means and 
  variances (\code{wt.moments}), or single number (\code{wt.var})
}


\author{
  Korbinian Strimmer (\url{https://strimmerlab.github.io}).
}

\seealso{\code{\link{weighted.mean}}, \code{\link{cov.wt}}.}

\examples{
# load corpcor library
library("corpcor")

# generate some data
p = 5
n = 5
X = matrix(rnorm(n*p), nrow = n, ncol = p)
w = c(1,1,1,3,3)/9


# standardize matrix
scale(X)
wt.scale(X)
wt.scale(X, w) # take into account data weights

}
\keyword{multivariate}
