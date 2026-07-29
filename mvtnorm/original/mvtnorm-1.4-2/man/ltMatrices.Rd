\name{ltMatrices}
\alias{ltMatrices}
\alias{syMatrices}
\alias{Tcrossprod}
\alias{Crossprod}
\alias{tcrossprod.ltMatrices}
\alias{crossprod.ltMatrices}
\alias{tcrossprod.syMatrices}
\alias{crossprod.syMatrices}
\alias{diagonals}
\alias{Mult}
\alias{Mult.ltMatrices}
\alias{Mult.syMatrices}
\alias{solve.ltMatrices}
\alias{logdet}
\alias{as.array.ltMatrices}
\alias{as.array.syMatrices}
\alias{diagonals.ltMatrices}
\alias{diagonals.syMatrices}
\alias{diagonals.matrix}
\alias{diagonals<-}
\alias{diagonals<-.ltMatrices}
\alias{diagonals<-.syMatrices}
\alias{diagonals.integer}
\alias{is.ltMatrices}
\alias{is.syMatrices}
\alias{as.ltMatrices}
\alias{as.ltMatrices.ltMatrices}
\alias{as.ltMatrices.syMatrices}
\alias{as.syMatrices}
\alias{is.chol}
\alias{is.invchol}
\alias{as.chol}
\alias{as.invchol}
\alias{Lower_tri}
\alias{chol.syMatrices}
\alias{aperm.chol}
\alias{aperm.invchol}
\alias{aperm.ltMatrices}
\alias{aperm.syMatrices}
\alias{deperma}
\alias{adddiag}
\alias{chol2cov}
\alias{invchol2chol}
\alias{chol2invchol}
\alias{invchol2cov}
\alias{invchol2pre}
\alias{chol2pre}
\alias{Dchol}
\alias{invcholD}
\alias{chol2cor}
\alias{invchol2cor}
\alias{invchol2pc}
\alias{chol2pc}
\alias{cov2chol}
\alias{cov2invchol}
\alias{invchol}
\alias{invchol.syMatrices}
\alias{vectrick}
\alias{standardize}
\alias{destandardize}
\title{
  Multiple Lower Triangular or Symmetric Matrices
}
\description{
  A class representing multiple lower triangular or symmetric matrices and some methods.
}
\usage{
ltMatrices(object, diag = FALSE, byrow = FALSE, names = TRUE)
syMatrices(object, diag = FALSE, byrow = FALSE, names = TRUE)
\S3method{as.array}{ltMatrices}(x, symmetric = FALSE, \dots)
\S3method{as.array}{syMatrices}(x, \dots)
\S3method{diagonals}{ltMatrices}(x, \dots)
\S3method{diagonals}{syMatrices}(x, \dots)
\S3method{diagonals}{matrix}(x, \dots)
\S3method{diagonals}{integer}(x, \dots)
diagonals(x) <- value
\S3method{diagonals}{ltMatrices}(x) <- value
\S3method{diagonals}{syMatrices}(x) <- value
\S3method{solve}{ltMatrices}(a, b, transpose = FALSE, \dots)
\S3method{chol}{syMatrices}(x, \dots)
\S3method{aperm}{chol}(a, perm, \dots)
\S3method{aperm}{invchol}(a, perm, \dots)
\S3method{aperm}{ltMatrices}(a, perm, \dots)
\S3method{aperm}{syMatrices}(a, perm, \dots)
deperma(chol = solve(invchol), permuted_chol = solve(permuted_invchol), 
        invchol, permuted_invchol, perm, score_schol)
\S3method{Mult}{ltMatrices}(x, y, transpose = FALSE, \dots)
\S3method{Mult}{syMatrices}(x, y, \dots)
Tcrossprod(x, diag_only = FALSE)
Crossprod(x, diag_only = FALSE)
\S3method{tcrossprod}{ltMatrices}(x, y = NULL, ...)
\S3method{tcrossprod}{syMatrices}(x, y = NULL, ...)
\S3method{crossprod}{ltMatrices}(x, y = NULL, ...)
\S3method{crossprod}{syMatrices}(x, y = NULL, ...)
logdet(x)
Lower_tri(x, diag = FALSE, byrow = attr(x, "byrow"))
is.ltMatrices(x)
is.syMatrices(x)
as.ltMatrices(x)
\S3method{as.ltMatrices}{ltMatrices}(x)
\S3method{as.ltMatrices}{syMatrices}(x)
as.syMatrices(x)
is.chol(x)
is.invchol(x)
as.chol(x)
as.invchol(x)
chol2cov(x)
invchol2chol(x)
chol2invchol(x)
invchol2cov(x)
invchol2pre(x)
chol2pre(x)
Dchol(x, D = 1 / sqrt(Tcrossprod(x, diag_only = TRUE)))
invcholD(x, D = sqrt(Tcrossprod(solve(x), diag_only = TRUE)))
chol2cor(x)
invchol2cor(x)
chol2pc(x)
invchol2pc(x)
cov2invchol(x)
cov2chol(x)
invchol(x, ...)
\S3method{invchol}{syMatrices}(x, ...)
vectrick(C, S, A, transpose = c(TRUE, TRUE))
standardize(chol, invchol)
destandardize(chol = solve(invchol), invchol, score_schol)
as.ltMatrices(x)
}
\arguments{
  \item{object}{a \code{matrix} representing the lower triagular elements of
                \eqn{N} lower triangular matrix, each of dimension \eqn{J \times
                J}. Dimensions of \code{object} depend on \code{diag}: With
                diagonal elements, \code{object} is a \eqn{J(J+1)/2 \times N} matrix,
                otherwise, the number of rows is \eqn{J(J - 1) / 2}.
}
  \item{diag}{logical, \code{object} contains diagonal elements if
              \code{TRUE}, otherwise unit diagonal elements are assumed.
}
  \item{byrow}{logical, \code{object} represents matrices in row-major order
               if \code{TRUE} or, otherwise, in column-major order.
}
  \item{names}{logical or character vector of length \eqn{J}.
}
  \item{symmetric}{logical, object is interpreted as a symmetric matrix if
\code{TRUE}.
}
  \item{diag_only}{logical, compute diagonal elements of crossproduct only
if \code{TRUE}.
}
  \item{x,chol,invchol,permuted_chol,permuted_invchol}{object of class \code{ltMatrices} or \code{syMatrices} (for \code{chol}).
}
  \item{value}{a matrix of diagonal elements to be assigned (of dimension \eqn{J \times N}).
}
  \item{a}{object of class \code{ltMatrices}.
}
  \item{perm}{a permutation of the covariance matrix corresponding to \code{a}.
}
  \item{D}{a matrix (of dimension \eqn{J \times N}) of diagonal elements to be multiplied with.
}
  \item{y}{matrix with \eqn{J} or \eqn{J^2} rows, interpreted as matrix of
           vectors or matrix of vectorised quadratic matrices.
}
  \item{b}{matrix with \eqn{J} rows.
}
  \item{C}{an object of class \code{ltMatrices}.}
  \item{S}{an object of class \code{ltMatrices} or a matrix with \eqn{J^2} rows 
           representing multiple \eqn{J x J} matrices (columns of vec operators).}
  \item{A}{an object of class \code{ltMatrices}.}
  \item{transpose}{a logical of length two indicating if \code{A} or \code{B} shall be transposed in \code{vectrick}.
                   For \code{solve}, this argument being true computes
                   \code{solve(t(a), b)} (in absence of a \code{t()} method for \code{ltMatrices}
                   objects).}
  \item{score_schol}{score matrix for a standardized \code{chol} object.}
  \item{\dots}{additional arguments, currently ignored.
}
}
\details{
\code{ltMatrices} interprets a matrix as lower triangular elements of
multiple lower triangular matrices. The corresponding class can be used to
store such matrices efficiently. Matrix multiplications, solutions to linear
systems, explicite inverses, and crossproducts can be computed based on such
objects. Details can be found in the \code{lmvnorm_src} package vignette.

\code{syMatrices} only store the lower triangular parts of multiple
symmetric matrices.
}
\value{
The constructor \code{ltMatrices} returns objects of class \code{ltMatrices}
with corresponding methods. The constructor \code{syMatrices} returns objects of class
\code{syMatrices} with a reduced set of methods.
}
\seealso{\code{vignette("lmvnorm_src", package = "mvtnorm")}}
\examples{

  J <- 4L
  N <- 2L
  dm <- paste0("d", 1:J)
  xm <- paste0("x", 1:N)
  (C <- ltMatrices(matrix(runif(N * J * (J + 1) / 2), 
                          ncol = N, dimnames = list(NULL, xm)), 
                   diag = TRUE, names = dm))

  ## dimensions and names
  dim(C)
  dimnames(C)
  names(C)

  ## subset
  C[,2:3]

  ## multiplication
  y <- matrix(runif(N * J), nrow = J)
  Mult(C, y)
  C %*% y

  ## solve
  solve(C)
  solve(C, y)

  ## tcrossprod
  Tcrossprod(C)
  tcrossprod(C)

  ## convert to matrix
  as.array(solve(C[1,]))[,,1]

}
\keyword{matrix}
