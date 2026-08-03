\name{pseudoinverse}
\alias{pseudoinverse}

\title{Pseudoinverse of a Matrix}
\usage{
pseudoinverse(m, tol)
}
\arguments{
  \item{m}{matrix}
  \item{tol}{tolerance - singular values larger than
             tol are considered non-zero (default value:
	     \code{tol = max(dim(m))*max(D)*.Machine$double.eps})
	     }
}
\description{
  The standard definition for the inverse of a matrix fails 
  if the matrix is not square or singular. However, one can
  generalize the inverse using singular value decomposition.
  Any rectangular real matrix M can be decomposed as

                  \deqn{M = U D V^{'},}{M = U D V',}

  where U and V are orthogonal, V' means V transposed, and 
  D is a diagonal matrix containing only the positive singular values
  (as determined by \code{tol}, see also \code{\link{fast.svd}}). 
  
   The pseudoinverse, also known as Moore-Penrose or generalized inverse
   is then obtained as
  
                  \deqn{iM = V D^{-1} U^{'}}{iM = V D^(-1) U' .}
			
}

\details{  
  The pseudoinverse has the property that the sum of the squares of all
  the entries in \code{iM \%*\% M - I}, where I is an appropriate
  identity matrix, is minimized. For non-singular matrices the
  pseudoinverse is equivalent to the standard inverse.
}

\value{
  A matrix (the pseudoinverse of m).
}


\author{
  Korbinian Strimmer (\url{https://strimmerlab.github.io}).
}
\seealso{\code{\link{solve}}, \code{\link{fast.svd}}}

\examples{
# load corpcor library
library("corpcor")

# a singular matrix
m = rbind(
c(1,2),
c(1,2)
)

# not possible to invert exactly
try(solve(m))

# pseudoinverse
p = pseudoinverse(m)
p

# characteristics of the pseudoinverse
zapsmall( m \%*\% p \%*\% m )  ==  zapsmall( m )
zapsmall( p \%*\% m \%*\% p )  ==  zapsmall( p )
zapsmall( p \%*\% m )  ==  zapsmall( t(p \%*\% m ) )
zapsmall( m \%*\% p )  ==  zapsmall( t(m \%*\% p ) )


# example with an invertable matrix
m2 = rbind(
c(1,1),
c(1,0)
)
zapsmall( solve(m2) ) == zapsmall( pseudoinverse(m2) )
}
\keyword{algebra}
