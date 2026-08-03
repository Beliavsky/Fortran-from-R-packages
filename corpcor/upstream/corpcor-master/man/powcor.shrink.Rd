\name{powcor.shrink}
\alias{powcor.shrink}
\alias{crossprod.powcor.shrink}

\title{Fast Computation of the Power of the Shrinkage Correlation Matrix}
\description{
  The function \code{powcor.shrink} efficiently computes the \code{alpha}-th power
  of the shrinkage correlation matrix produced by \code{\link{cor.shrink}}. 

  For instance, this function may be used for fast computation of the (inverse)
  square root of the  shrinkage correlation matrix (needed, e.g., for decorrelation).

  \code{crossprod.powcor.shrink} efficiently computes \eqn{R^{\alpha} y} without
  actually computing the full matrix \eqn{R^{\alpha}}.
}
\usage{
powcor.shrink(x, alpha, lambda, w, verbose=TRUE)
crossprod.powcor.shrink(x, y, alpha, lambda, w, verbose=TRUE)
}
\arguments{
  \item{x}{a data matrix}
  \item{y}{a matrix, the number of rows of y must be the same as the number of columns of x}
  \item{alpha}{exponent}
  \item{lambda}{the correlation shrinkage intensity (range 0-1).
                If \code{lambda} is not specified (the default) it is estimated
		using an analytic formula from Sch\"afer and Strimmer (2005) 
		- see \code{\link{cor.shrink}}.  
		For \code{lambda=0} the empirical correlations are recovered.}
  \item{w}{optional: weights for each data point - if not specified uniform weights are assumed
            (\code{w = rep(1/n, n)} with \code{n = nrow(x)}).}	
  \item{verbose}{output status while computing (default: TRUE)}
}

\details{
  This function employs a special matrix identity to speed up the computation of
  the matrix power of the shrinkage correlation matrix (see Zuber and Strimmer 2009 for details).

  Apart from a scaling factor the shrinkage correlation matrix computed 
  by \code{\link{cor.shrink}}  takes
  on the form

   \deqn{Z = I_p + V M V^T ,}   

  where \code{V M V^T} is a multiple of the empirical correlation matrix.
  Crucially, \code{Z} is a matrix of size \code{p} times \code{p}
   whereas  \code{M} is a potentially much smaller matrix of size \code{m} times \code{m},
   where \code{m} is the rank of the empirical correlation matrix.  

  In order to calculate the \code{alpha}-th power of \code{Z}
  the function uses the identity

    \deqn{Z^\alpha = I_p - V (I_m -(I_m + M)^\alpha) V^T}

  requiring only the computation of the \code{alpha}-th power of the \code{m} by \code{m} matrix
  \eqn{I_m + M}.  This trick enables substantial computational savings especially when the number
  of observations is much smaller than the number of variables.
   
   Note that the above identity is related but not identical to the Woodbury matrix
   identity for inversion of a matrix.
   For \eqn{\alpha=-1} the above identity reduces to

    \deqn{Z^{-1} = I_p - V (I_m -(I_m + M)^{-1}) V^T ,}

   whereas the Woodbury matrix identity equals

    \deqn{Z^{-1} = I_p - V (I_m + M^{-1})^{-1} V^T .}

}

\value{
  \code{powcor.shrink} returns a matrix of the same size as the correlation matrix \code{R}

 \code{crossprod.powcor.shrink} returns a matrix of the same size as \code{R} \code{y}.
    
}
\author{
  Verena Zuber, A. Pedro Duarte Silva, and Korbinian Strimmer (\url{https://strimmerlab.github.io}).
  
}

\references{

  Zuber, V., and K. Strimmer. 2009. Gene ranking and biomarker
  discovery under correlation.  Bioinformatics \bold{25}:2700-2707.
  <DOI:10.1093/bioinformatics/btp460>

  Zuber, V.,  A. P. Duarte Silva, and K. Strimmer. 2012. A novel algorithm for 
  simultaneous SNP selection in high-dimensional genome-wide association studies.
  BMC Bioinformatics 13: 284
  <DOI:10.1186/1471-2105-13-284>

}


\seealso{\code{\link{invcor.shrink}}, \code{\link{cor.shrink}}, \code{\link{mpower}}.}

\examples{
# load corpcor library
library("corpcor")

# generate data matrix
p = 500
n = 10
X = matrix(rnorm(n*p), nrow = n, ncol = p)

lambda = 0.23  # some arbitrary lambda

### computing the inverse ###
# slow
system.time(
  (W1 = solve(cor.shrink(X, lambda=lambda)))
)

# very fast
system.time(
  (W2 = powcor.shrink(X, alpha=-1, lambda=lambda))
)

# no difference
sum((W1-W2)^2)

### computing the square root ###

system.time(
  (W1 = mpower(cor.shrink(X, lambda=lambda), alpha=0.5))
)

# very fast
system.time(
  (W2 = powcor.shrink(X, alpha=0.5, lambda=lambda))
)

# no difference
sum((W1-W2)^2)


### computing an arbitrary power (alpha=1.23) ###

system.time(
  (W1 = mpower(cor.shrink(X, lambda=lambda), alpha=1.23))
)

# very fast
system.time(
  (W2 = powcor.shrink(X, alpha=1.23, lambda=lambda))
)

# no difference
sum((W1-W2)^2)


### fast computation of cross product

y = rnorm(p)

system.time(
  (CP1 = crossprod(powcor.shrink(X, alpha=1.23, lambda=lambda), y))
)

system.time(
  (CP2 = crossprod.powcor.shrink(X, y, alpha=1.23, lambda=lambda))
)

# no difference
sum((CP1-CP2)^2)

}
\keyword{multivariate}
