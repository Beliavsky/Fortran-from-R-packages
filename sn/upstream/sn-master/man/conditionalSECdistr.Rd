%  file sn/man/conditionalSECdistr.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{conditionalSECdistr}
\alias{conditionalSECdistr}

\title{Skew-normal conditional distribution}
 
\description{For a multivariate (extended) skew-normal distribution, compute 
its conditional distribution for given values of some of its components.}

\usage{
  conditionalSECdistr(object, fixed.comp, fixed.values, name, drop = TRUE)
}

\arguments{

  \item{object}{an object of class \code{SECdistrMv} with \code{family="SN"}
     or \code{family="ESN"}. }
 
  \item{fixed.comp}{a vector containing a subset of \code{1:d} which selects
    the components whose values are to be fixed, if \code{d} denotes the
    dimensionality of the distribution.}

  \item{fixed.values}{a numeric vector of values taken on by the components
    \code{fixed.comp}; it must be of the same length of \code{fixed.comp}.}
 
  \item{name}{an optional character string with the name of the outcome 
    distribution; if missing, one such string is constructed.}

  \item{drop}{logical (default=\code{TRUE}), to indicate whether the
    returned object must be of class \code{SECdistrUv} when
    \code{length(fixed.comp)+1=d}.}

}

\value{an object of class \code{SECdistrMv}, except in the case when 
  \code{drop=TRUE} operates, leading to an object of class 
  \code{SECdistrUv-class}.}
 
\details{For background information, see Section 5.3.2 of the reference below.}

\references{
  Azzalini, A. and Capitanio, A. (2014). \emph{The Skew-normal and Related 
  Families}. Cambridge University Press, IMS Monographs series.
}

\seealso{\code{\link{makeSECdistr}}, \code{\link{SECdistrMv-class}}, 
  \code{\link{affineTransSECdistr}} }

\examples{
Omega <- diag(3) + outer(1:3,1:3)
sn <- makeSECdistr(dp=list(xi=rep(0,3), Omega=Omega, alpha=1:3), family="SN")
esn <- conditionalSECdistr(sn, fixed.comp=2, fixed.values=1.5)
show(esn)
}
\keyword{multivariate}
\keyword{distribution}
