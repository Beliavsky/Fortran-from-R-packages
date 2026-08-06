%  file sn/man/extractSECdistr.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{extractSECdistr}
\alias{extractSECdistr}
\concept{skew-elliptical distribution}
\title{Extract the SEC error distribution from an object created by \code{selm}}

\description{
  Given an object created by a call to \code{selm}, the function delivers
  the \acronym{SEC} distribution representing the stochastic term of the 
  fitted model
}

\usage{
  extractSECdistr(object, name, compNames)
}

\arguments{
  \item{object}{an object of class \code{selm} or \code{mselm}, 
    as created by \code{\link{selm}}.}
 
  \item{name}{an optional character string representing the name of the 
    outcome distribution; if missing, a string is constructed from the
    \code{object} ingredients.}

  \item{compNames}{in the multivariate case, an optional vector  of character 
    strings  with the names of the components of the error distribution; 
    if missing, one such vector is constructed from the \code{object}
    ingredients.}      
}

\value{An object of class \code{SECdistrMv} or \code{SECdistrUv}, 
  depending of the class of \code{object}.}
 
\section{Details}{
When the formula of the fitted model includes only the constant \code{1},
the returned object represents the fitted \acronym{SEC} distribution.
If the formula includes additional terms, the linear predictor is eliminated
and the returned object corresponds to the error term of the model; 
hence the location parameter \code{xi} in the \acronym{DP} parameterization
is set to zero.

The returned object can be submitted to tools
available for objects created by \code{\link{makeSECdistr}}, 
such as \code{\link{summary.SECdistr}}, \code{\link{conditionalSECdistr}} and
and so on.}
 
\seealso{\code{\link{selm}}, \code{\link{makeSECdistr}}}
  
\examples{
data(ais)
m2 <- selm(log(Fe) ~ 1, family="ST", data=ais, fixed=list(nu=8))
f2 <- extractSECdistr(m2)
show(f2)
#
m4 <- selm(cbind(BMI, LBM) ~ 1, family="SN", data=ais)
f4 <- extractSECdistr(m4)
mean(f4)
vcov(f4)
}

\keyword{multivariate}
\keyword{distribution}
