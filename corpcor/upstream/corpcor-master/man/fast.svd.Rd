\name{fast.svd}
\alias{fast.svd}

\title{Fast Singular Value Decomposition}

\description{
  \code{fast.svd} returns the singular value decomposition of
  a rectangular real matrix 
  
                             \deqn{M = U D V^{'},}{M = U D V',}

  where \eqn{U} and \eqn{V} are orthogonal matrices with \eqn{U' U = I}
  and \eqn{V' V = I}, and \eqn{D} is a diagonal matrix containing the 
  singular values (see \code{\link{svd}}).
  
  The main difference to the native version \code{\link{svd}} is that 
  \code{fast.svd} is substantially faster for  "fat" (small n, large p)
  and "thin" (large n, small p) matrices.
  In this case the decomposition of \eqn{M} can be
  greatly sped up by first computing the SVD of either \eqn{M M'} (fat matrices) or
  \eqn{M' M} (thin matrices), rather than that of \eqn{M}.
  
  A second difference to \code{\link{svd}} is that \code{fast.svd} only
  returns the \emph{positive} singular values (thus the dimension of \eqn{D}
  always equals the rank of \eqn{M}).  Note that the singular
  vectors computed by \code{fast.svd} may differ in sign from those computed 
  by \code{\link{svd}}. 
}



\usage{
fast.svd(m, tol)
}
\arguments{
  \item{m}{matrix}
  \item{tol}{tolerance - singular values larger than
             tol are considered non-zero (default value:
	     \code{tol = max(dim(m))*max(D)*.Machine$double.eps})
	     }
}

\details{
  For "fat" \eqn{M} (small n, large p) the SVD decomposition of \eqn{M M'} yields
  
            \deqn{M M^{'} = U D^2 U}{M M' = U D^2 U'}
	     
  As the matrix \eqn{M M'} has dimension n x n only, this is faster to compute
  than SVD of \eqn{M}.  The \eqn{V} matrix is subsequently obtained by
  
             \deqn{V = M^{'} U D^{-1}}{V = M' U D^(-1)}  
	            
  Similarly, for "thin" \eqn{M} (large n, small p), the decomposition of \eqn{M' M}
  yields
  
             \deqn{M^{'} M = V D^2 V^{'}}{M' M = V D^2 V'}
	    
  which is also quick to compute as \eqn{M' M} has only dimension p x p.  The 
  \eqn{U} matrix is then computed via
  
             \deqn{U = M V D^{-1}}{U = M V D^(-1)}
}

\value{
  A list with the following components:
  \item{d}{a vector containing the \emph{positive} singular values}
  \item{u}{a matrix with the corresponding left singular vectors}
  \item{v}{a matrix with the corresponding right singular vectors}
}


\author{
  Korbinian Strimmer (\url{https://strimmerlab.github.io}).
}
\seealso{\code{\link{svd}}, \code{\link{solve}}.}

\examples{
# load corpcor library
library("corpcor")


# generate a "fat" data matrix
n = 50
p = 5000
X = matrix(rnorm(n*p), n, p)

# compute SVD
system.time( (s1 = svd(X)) ) 
system.time( (s2 = fast.svd(X)) )


eps = 1e-10
sum(abs(s1$d-s2$d) > eps)
sum(abs(abs(s1$u)-abs(s2$u)) > eps)
sum(abs(abs(s1$v)-abs(s2$v)) > eps)
}
\keyword{algebra}
