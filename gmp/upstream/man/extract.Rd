\name{extract}
\title{Extract or Replace Parts of a 'bigz' or 'bigq' Object}
\alias{[.bigz}
%F\alias{[<-}% otherwise warning " Objects in \usage w/o \alias : '[<-' "
\alias{[<-.bigz}
\alias{[[.bigz}
\alias{[[<-.bigz}
\alias{c.bigz}
\alias{rep.bigz}
\alias{length.bigz}
\alias{length<-.bigz}
\alias{[.bigq}
\alias{[<-.bigq}
\alias{[[.bigq}
\alias{[[<-.bigq}
%-
\alias{c.bigq}
\alias{rep.bigq}
\alias{length.bigq}
\alias{length<-.bigq}
\description{
  Operators acting on vectors, arrays and lists to extract or replace subsets.
}
\usage{
\method{[}{bigz}(x, i=NULL, j=NULL, drop = TRUE)
\method{[}{bigq}(x, i=NULL, j=NULL, drop = TRUE)
%F \method{[<-}{bigq}(x, i=NULL, j=NULL, value) %    ==> "the same", manually -- R check-bug:
%F x[i=NULL, j=NULL] <- value % still WARNING .. usage in doc.. .. but not in code: `[<-`
##___ In the following, only the bigq method is mentioned (but 'bigz' is "the same"): ___
\method{c}{bigq}(\dots, recursive = FALSE)
\method{rep}{bigq}(x, times=1, length.out=NA, each=1, \dots)
}
\arguments{
  \item{x}{\R object of class \code{"bigz"} or \code{"bigq"}, respectively.}
  \item{\dots}{further arguments, notably for \code{c()}.}
  \item{i,j}{indices, see standard \R subsetting and subassignment.}
  \item{drop}{logical, unused here, i.e., matrix subsetting \bold{always} returns a matrix, here!}
%F \item{value}{\R object, typically of same \code{\link{class}} as
%F   \code{x}, or also \code{\link{numeric}}.}
  \item{times, length.out, each}{integer; typically only \emph{one} is
    specified; for more see \code{\link{rep}} (standard \R, package \pkg{base}).
  }
  \item{recursive}{from \code{c()}'s default method; disregarded here}
}

\examples{
  a <- as.bigz(123)
  ## indexing "outside" --> extends the vectors (filling with NA)
  a[2] <- a[1]
  a[4] <- -4

  ## create a vector of 3 a
  c(a,a,a)

  ## repeate a 5 times
  rep(a,5)

  ## with matrix: 3 x 2
  m <- matrix.bigz(1:6,3)

  m[1,] # the first row
  m[1,, drop=TRUE] # the same: drop does *not* drop
  m[1]
  m[-c(2,3),]
  m[-c(2,3)]
  m[c(TRUE,FALSE,FALSE)]
\dontshow{
  stopifnot(identical(c(a,a), rep(a,2)),
            identical(rep(a,0), a[FALSE]), length(a[FALSE]) == 0,
            identical(m[1,], m[-(2:3),]),
            identical(m[1,], m[c(TRUE,FALSE,FALSE),])
  )
}
  ##_modification on matrix
  m[2,-1] <- 11
}
\keyword{arith}
