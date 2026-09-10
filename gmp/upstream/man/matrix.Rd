\name{matrix}
\alias{matrix}
\alias{matrix.default}
\alias{matrix.bigz}
\alias{matrix.bigq}
\alias{is.matrixZQ}
\alias{as.matrix.bigz}
\alias{as.matrix.bigq}
\alias{as.vector.bigq}
\alias{as.vector.bigz}
%-- if this worked with S4:
% \alias{\%*\%,bigq,bigq-method}
% \alias{\%*\%,bigq,ANY-method}
% \alias{\%*\%,ANY,bigq-method}
% \alias{\%*\%,bigz,bigz-method}
% \alias{\%*\%,bigz,ANY-method}
% \alias{\%*\%,ANY,bigz-method}
% \alias{crossprod,bigq,bigq-method}
% \alias{crossprod,bigq,ANY-method}
% \alias{crossprod,ANY,bigq-method}
% \alias{crossprod,bigz,bigz-method}
% \alias{crossprod,bigz,ANY-method}
% \alias{crossprod,ANY,bigz-method}
% \alias{tcrossprod,bigq,bigq-method}
% \alias{tcrossprod,bigq,ANY-method}
% \alias{tcrossprod,ANY,bigq-method}
% \alias{tcrossprod,bigz,bigz-method}
% \alias{tcrossprod,bigz,ANY-method}
% \alias{tcrossprod,ANY,bigz-method}
%-- need to stay S3 -- as primitive %*% does not dispatch on setOldClass()ed pseudo-S4
\alias{\%*\%}
\alias{\%*\%.default}
\alias{\%*\%.bigq}
\alias{\%*\%.bigz}
\alias{crossprod}
\alias{crossprod.default}
\alias{crossprod.bigq}
\alias{crossprod.bigz}
\alias{tcrossprod}
\alias{tcrossprod.default}
\alias{tcrossprod.bigq}
\alias{tcrossprod.bigz}
%
\alias{ncol.bigq}
\alias{ncol.bigz}
\alias{nrow.bigq}
\alias{nrow.bigz}
\alias{cbind.bigz}
\alias{cbind.bigq}
\alias{rbind.bigz}
\alias{rbind.bigq}
\alias{t.bigq}
\alias{t.bigz}
\alias{dim.bigq}
\alias{dim<-.bigq}
\alias{dim.bigz}
\alias{dim<-.bigz}

\title{Matrix manipulation with gmp}
\description{
  Overload of \dQuote{all} standard tools useful for matrix manipulation adapted
  to large numbers.
}
% \S4method{\%*\%}{bigz,bigz}(x, y)
% \S4method{\%*\%}{bigz,ANY}(x, y)
% \S4method{\%*\%}{ANY,bigz}(x, y)
\usage{
\S3method{matrix}{bigz}(data = NA, nrow = 1, ncol = 1, byrow = FALSE, dimnames = NULL, mod = NA,...)

is.matrixZQ(x)

\S3method{\%*\%}{bigz}(x, y)
\S3method{\%*\%}{bigq}(x, y)
\S3method{crossprod}{bigq}(x, y=NULL,...)
\S3method{tcrossprod}{bigz}(x, y=NULL,...)

\S3method{cbind}{bigz}(..., deparse.level=1)
\S3method{rbind}{bigq}(..., deparse.level=1)
## ..... etc
}
\arguments{
  \item{data}{an optional data vector}
  \item{nrow}{the desired number of rows }
  \item{ncol}{the desired number of columns}
  \item{byrow}{logical.  If \code{FALSE} (the default), the matrix is filled by
    columns, otherwise the matrix is filled by rows.}
  \item{dimnames}{not implemented for \code{"bigz"} or \code{"bigq"} matrices.}
  \item{mod}{optional modulus (when \code{data} is \code{"bigz"}).}

  \item{x,y}{numeric, \code{bigz}, or \code{bigq} matrices or vectors.}

  \item{..., deparse.level}{arguments from the generic; \emph{not} made use
  of, i.e., disregarded here.}
}
\details{
  The extract function (\code{"["}) is the same use for vector or
  matrix.  Hence, \code{x[i]} returns the same values as \code{x[i,]}.
  This is not considered a feature and may be changed in the future
  (with warnings).

  All matrix multiplications should work as with numeric matrices.

  Special features concerning the \code{"\link{bigz}"} class: the
  modulus can be
  \describe{
    \item{Unset:}{Just play with large numbers}
    \item{Set with a vector of size 1:}{Example:
      \code{matrix.bigz(1:6,nrow=2,ncol=3,mod=7)}
      This means you work
      in \eqn{Z/nZ}, for the whole matrix.  It is the only case
      where the \code{\%*\%} and \code{\link{solve}} functions will work
      in \eqn{Z/nZ}.}
    \item{Set with a vector smaller than data:}{Example:
      \code{matrix.bigz(1:6,nrow=2,ncol=3,mod=1:5)}. Then, the modulus
      is repeated to the end of data.  This can be used to define a
      matrix with a different modulus at each row.
    }
    \item{Set with same size as data:}{Modulus is defined for each cell}
  }
}
\value{
  \code{matrix()}: A matrix of class \code{"bigz"} or \code{"bigq"}.

  \code{is.matrixZQ()}: \code{\link{TRUE}} or \code{FALSE}.

  \code{dim()}, \code{ncol()}, etc: integer or \code{NULL}, as for
  simple matrices.

  \code{cbind(x,y,...)} and \code{rbind(x,y,...)} now (2024-01, since
  \pkg{gmp} version 0.9-5), do drop \code{deparse.level=.} instead of
  wrongly creating an extra column or row \emph{and} the \code{"bigz"}
  method takes all arguments into account and calls the \code{"bigq"}
  method in case of arguments inheriting from \code{"\link{bigq}"}.
}
\author{Antoine Lucas and Martin Maechler}

\seealso{Solving a linear system: \code{\link{solve.bigz}}.
  \code{\link[base]{matrix}}
}

\examples{
V <- as.bigz(v <- 3:7)
crossprod(V)# scalar product
(C <- t(V))
stopifnot(dim(C) == dim(t(v)), C == v,
          dim(t(C)) == c(length(v), 1),
          crossprod(V) == sum(V * V),
         tcrossprod(V) == outer(v,v),
          identical(C, t(t(C))),
          is.matrixZQ(C), !is.matrixZQ(V), !is.matrixZQ(5)
	)

## a matrix
x <- diag(1:4)
## invert this matrix
(xI <- solve(x))

## matrix in Z/7Z
y <- as.bigz(x,7)
## invert this matrix (result is *different* from solve(x)):
(yI <- solve(y))
stopifnot(yI \%*\% y == diag(4),
          y \%*\% yI == diag(4))

## matrix in Q
z  <- as.bigq(x)
## invert this matrix (result is the same as solve(x))
(zI <- solve(z))

stopifnot(abs(zI - xI) <= 1e-13,
          z \%*\% zI == diag(4),
          identical(crossprod(zI), zI \%*\% t(zI))
         )

A <- matrix(2^as.bigz(1:12), 3,4)
for(a in list(A, as.bigq(A, 16), factorialZ(20), as.bigq(2:9, 3:4))) {
  a.a <- crossprod(a)
  aa. <- tcrossprod(a)
  stopifnot(identical(a.a, crossprod(a,a)),
 	    identical(a.a, t(a) \%*\% a)
            ,
            identical(aa., tcrossprod(a,a)),
	    identical(aa., a \%*\% t(a))
 	   )
}# {for}
}
\keyword{arith}
