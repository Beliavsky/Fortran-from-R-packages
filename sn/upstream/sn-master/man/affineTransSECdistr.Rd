%  file sn/man/affineTransSECdistr.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{affineTransSECdistr}
\encoding{UTF-8}

\alias{marginalSECdistr}
\alias{affineTransSECdistr}
\title{Affine transformations and marginals of a skew-elliptical distribution}

\description{
  Given a multivariate random variable \eqn{Y} with skew-elliptical
 (\acronym{SEC}) distribution, compute the distribution 
  of a (possibly multivariate) marginal or the distribution
  of an affine  transformation \eqn{a + A^{\top}Y}{a + A'Y}.}

\usage{
  affineTransSECdistr(object, a, A, name, compNames, drop=TRUE)
  marginalSECdistr(object, comp, name, drop=TRUE)
}

\arguments{
  \item{object}{an object of class \code{SECdistrMv} which identifies the 
    source random variable, as created by
    \code{\link{makeSECdistr}} or by \code{\link{extractSECdistr}} or 
    by a previous call to these functions}

  \item{a}{a numeric vector with the length \code{ncol(A)}.}
  
  \item{A}{a full-rank matrix with \code{nrow(A)} equal to the dimensionality
    \kbd{d} of the random variable identified by \code{object}. }
    
  \item{name}{an optional character string representing the name of the 
    outcome distribution; if missing, one such string is constructed.}

  \item{compNames}{an optional vector of length \code{ncol(A)} of character 
    strings  with the names of the components of the outcome distribution; 
    if missing, one such vector is constructed.}
   
  \item{drop}{a logical flag (default value: \code{TRUE}), operating only if
    the returned object has dimension \code{d=1}, in which case it indicates
     whether this object must be of class \code{SECdistrUv}.}
  
  \item{comp}{a vector formed by a subset of \code{1:d} which indicates which
    components must be extracted from \code{object}, on denoting by \code{d} 
    its dimensionality.}
      
}
\value{If \code{object} defines the distribution of a \acronym{SEC} random
 variable \eqn{Y}, \code{affineTransSECdistr} computes  the 
 distribution of \eqn{a+A'Y} and \code{marginalSECdistr} computes the marginal
 distribution of the \code{comp} components. In both cases the returned
 object is of class \code{SECdistrMv}, except when \code{drop=TRUE}
 operates, leading to an object of class \code{SECdistrUv}.}

\section{Background}{These functions implement formulae given in Sections 
5.1.4, 5.1.6 and 6.2.2 of the reference below.}

\references{
 Azzalini, A. with the collaboration of Capitanio, A. (2014). 
 \emph{The Skew-Normal and Related Families}. 
 Cambridge University Press, IMS Monographs series. 
}  
 
\seealso{\code{\link{makeSECdistr}}, \code{\link{extractSECdistr}},
  \code{\link{SECdistrMv-class}}}
  
\examples{
dp3 <- list(xi=1:3, Omega=toeplitz(1/(1:3)), alpha=c(3,-1,2), nu=5)
st3 <- makeSECdistr(dp3, family="ST", name="ST3", compNames=c("U", "V", "W"))
A <- matrix(c(1,-1,1, 3,0,-2), 3, 2)
new.st <- affineTransSECdistr(st3, a=c(-3,0), A=A)
#
st2 <- marginalSECdistr(st3, comp=c(3,1), name="2D marginal of ST3")
}

\keyword{multivariate}
\keyword{distribution}
