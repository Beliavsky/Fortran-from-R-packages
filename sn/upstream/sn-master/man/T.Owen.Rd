%  file sn/man/T.Owen.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 1997-2013 Adelchi Azzalini
%---------------------
\name{T.Owen}
\alias{T.Owen}
\title{
Owen's function
}
\description{Evaluates function \eqn{T(h,a)} studied by D.B.Owen}

\usage{
T.Owen(h, a, jmax=50, cut.point=8)
}

\arguments{

\item{h}{a numeric vector. Missing values (\code{NA}s) and \code{Inf} are
 allowed.}

\item{a}{a numeric value. \code{Inf} is allowed.}

\item{jmax}{an integer scalar value which regulates the accuracy of the result.
  See Section \sQuote{Details} below for explanation.
}
\item{cut.point}{a scalar value which regulates the behaviour of the algorithm,
 as explained in Section \sQuote{Details} below (default value: \code{8}).}
 }

\value{a numeric vector.}

\details{
If \code{a>1} and \code{0<h<=cut.point}, a series expansion is used,
truncated after \code{jmax} terms.
If \code{a>1} and \code{h>cut.point}, an asymptotic approximation is used.
In the other cases, various reflection properties of the function
are exploited. See the reference below for more information.
}

\section{Background}{
The function \emph{T(h,a)} studied by Owen (1956) is useful for the computation 
of the bivariate normal distribution function and related quantities, 
including the distribution function of a skew-normal variate; see \code{psn}.
See the reference below for more information on function \eqn{T(h,a)}.
}

\author{Adelchi Azzalini and Francesca Furlan}

\references{
Owen, D. B. (1956).
Tables for computing bivariate normal probabilities.
\emph{Ann. Math. Statist.} \bold{27}, 1075-1090.
}

\seealso{ \code{\link{psn}}}

\examples{ owen <- T.Owen(1:10, 2)}

\keyword{math}
