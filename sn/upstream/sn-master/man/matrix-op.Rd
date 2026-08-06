%  file sn/man/vech.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 1998,2013 Adelchi Azzalini
%---------------------
\name{matrix-op}
\encoding{UTF-8}
\alias{vech}
\alias{vech2mat}
\alias{duplicationMatrix}
\alias{tr}
\alias{blockDiag}

\title{vech, tr and other matrix operators}

\description{vech and other matrix operators}
  
\usage{
vech(A)
vech2mat(v)
duplicationMatrix(n)
tr(x)
blockDiag(...)
}

\arguments{
\item{A}{a (symmetric) square numeric matrix.}

\item{v}{a numeric vector such that \code{length(v)=n*(n+1)/2} for some
positive integer \code{n}.}

\item{n}{a positive integer number; default is \code{n=1}.}

\item{x}{a square numeric matrix.}

\item{...}{an abitrary numer of matrices or objects coercible into matrices.}
}

\value{a vector in case of \code{vech}, a scalar in case of \code{tr}, 
  otherwise a matrix.}

\section{Details}{%
For a  square matrix \code{A}, \code{vech(A)} returns the vector formed
by the lower triangular portion of the matrix, including the diagonal; 
usually, this only makes sense for a symmetric matrix of numeric values. 
If \code{v=vech(M)} where \code{M} is a symmetric numeric matrix, 
\code{vech2mat(v)} performs the inverse operation and returns the original 
matrix \code{M}; this explain the requirement on \code{length(v)}.

For a positive integer \code{n}, \code{D=duplicationMatrix(n)} is a matrix 
of dimension \code{(n^2, n*(n+1)/2)} such that \code{D \%*\% vech(M)} returns 
the \code{vec}-form of a symmetric matrix \code{M} of 
order \code{n}, that is, the vector which stacks the columns of \code{M};
for more information, see Section 3.8 of Magnus and Neudecker (1988). 

For a square numeric matrix \code{x}, \code{tr(x)} returns its trace.

\code{blockDiag(...)} creates a block-diagonal matrix from a set of matrices 
or objects coercible into matrices. Generally, this is useful only for
numeric objects.
}

\section{Author}{Adelchi Azzalini; 
the original Octave code of \code{duplicationMatrix} is by Kurt Hornik.}
 
\references{
Magnus, Jan R. and Neudecker, Heinz (1988). \emph{Matrix differential
calculus with application in statistics and econometrics}. 
Wiley series in probability and statistics.
}

\examples{
M <- toeplitz(1:4)
v <- vech(M)
vech2mat(v) - M
D <- duplicationMatrix(ncol(M))
# D %*% vech(M) - as.vector(M), must be a one-column matrix of 0s
tr(outer(1:4,2:5))
blockDiag(M[1:2,], 1:2, diag(5:6))
} 

\keyword{math}
\concept{matrix operator}
\concept{duplication matrix}
