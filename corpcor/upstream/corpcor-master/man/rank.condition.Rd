\name{rank.condition}
\alias{is.positive.definite}
\alias{make.positive.definite}
\alias{rank.condition}

\title{Positive Definiteness of a Matrix, Rank and Condition Number}
\usage{
is.positive.definite(m, tol, method=c("eigen", "chol"))
make.positive.definite(m, tol)
rank.condition(m, tol)

}
\arguments{
  \item{m}{a matrix (assumed to be real and symmetric)}
  
  \item{tol}{tolerance for singular values and for absolute eigenvalues
              -  only those with values larger than
             tol are considered non-zero (default:
	     \code{tol = max(dim(m))*max(D)*.Machine$double.eps})
	     }
   \item{method}{Determines the method to check for positive definiteness:
     eigenvalue computation (\code{eigen}, default) or Cholesky decomposition
     (\code{chol}).}
     
}
\description{
  \code{is.positive.definite} tests whether all eigenvalues of a symmetric matrix
  are positive.  
  
  \code{make.positive.definite} computes the nearest positive definite of a
  real symmetric matrix, using the algorithm of NJ Higham (1988) <DOI:10.1016/0024-3795(88)90223-6>.

   \code{rank.condition} estimates the rank and the condition
  of a matrix by 
  computing its singular values D[i] (using  \code{\link{svd}}).
  The rank of the matrix is the number of singular values \code{D[i] > tol})
  and the condition is the ratio  of the largest and the smallest
  singular value.   
  
  The definition \code{tol= max(dim(m))*max(D)*.Machine$double.eps} 
  is exactly compatible with the conventions used in "Octave" or "Matlab".

  Also note that it is not checked whether the input matrix m is real and symmetric.

}

\value{
  \code{is.positive.definite} returns
  a logical value (\code{TRUE} or \code{FALSE}).
  
  
  \code{rank.condition} returns a list object with the following components:
  \item{rank}{Rank of the matrix.}
  \item{condition}{Condition number.}
  \item{tol}{Tolerance.}
  
  \code{make.positive.definite} returns a symmetric positive definite matrix.
}


\author{
  Korbinian Strimmer (\url{https://strimmerlab.github.io}).
}
\seealso{\code{\link{svd}}, \code{\link{pseudoinverse}}.}

\examples{
# load corpcor library
library("corpcor")

# Hilbert matrix
hilbert = function(n) { i = 1:n; 1 / outer(i - 1, i, "+") }

# positive definite ?
m = hilbert(8)
is.positive.definite(m)

# numerically ill-conditioned
m = hilbert(15)
rank.condition(m)

# make positive definite
m2 = make.positive.definite(m)
is.positive.definite(m2)
rank.condition(m2)
m2 - m
}
\keyword{algebra}
