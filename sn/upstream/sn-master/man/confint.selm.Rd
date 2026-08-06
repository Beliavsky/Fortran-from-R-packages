%  file sn/man/profile.selm.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2016 Adelchi Azzalini
%---------------------

\name{confint.selm}
\alias{confint.selm}
\alias{confint,selm-method}           
% \alias{confint.selm-class}
\concept{confidence interval}
\title{Confidence intervals for parameters of a selm-class object}

\description{
Computes confidence intervals for parameters in a selm-class object produces
by \code{selm} fit when the response variable is univariate.}
\usage{
  \S4method{confint}{selm}(object, parm, level=0.95, param.type, tol=1e-3, ...)

}

\arguments{
  \item{object}{an object of class \code{selm} as produced by a call to
    function \code{selm} with univariate response.}
  
  \item{parm}{a specification of which parameters are to be given confidence 
    intervals, either a vector of numbers or a vector of names. If missing, all 
    parameters are considered.}
 
  \item{level}{the confidence level required (default value is \code{0.95}).} 
  
  \item{param.type}{
    a character string with the required parameterization; it must be either
    \code{"CP"} or \code{"DP"} or \code{"pseudo-CP"}, or possibly their 
    equivalent lowercase.}
  
  \item{tol}{the desired accuracy (convergence tolerance); this is a parameter 
    passed to \code{\link[stats]{uniroot}} for computing the roots of the
    likelihood-based confidence interval for \code{alpha}.}
   
  \item{...}{not used, only there for compatibility reasons.}
  
}

\details{
A description of the methodology underlying \code{confint.selm} is provided 
in the technical note by Azzalini (2016).
This document also explains why in certain cases an interval is not 
produced and \code{NA}'s are returned as endpoint.}

\value{An invisible list whose components, described below, are partly 
different in the one- and the two-parameter cases.
  \item{call}{the calling statement}
  \item{<param1>}{values of the first parameter}
  \item{<param2>}{values of the second parameter (in a two-parameter case)}
  \item{logLik}{numeric vector or matrix of the profile log-likelihood values}
  \item{confint}{in the one-parameter case, the confidence interval}
  \item{level}{in the one-parameter case, the confidence level}
  \item{deviance.contour}{in the two-parameter case, a list of lists whose
    elements identify each curve of the contour plot}
  }
  

\author{Adelchi Azzalini}

\references{
 Azzalini, A. (2016). 
 Derivation of various types of intervals from a \code{selm} object.
 Technical note distributed with the documentation of the \R package
 \code{sn} in file 
 \href{../doc/selm-intervals.pdf}{\code{selm-intervals.pdf}}  within 
 section \sQuote{User guide, package vignettes and other documentation}.
 
% Azzalini, A. with the collaboration of Capitanio, A. (2014). 
% \emph{The Skew-Normal and Related Families}. 
% Cambridge University Press, IMS Monographs series. 
}
 
% \note{}
 

\seealso{
\code{\link{selm}}, \code{\link{summary.selm}}, \code{\link{profile.selm}},

\code{\link{makeSECdistr}} for the \acronym{CP}/\acronym{DP} parameterizations,

\code{\link[stats]{uniroot}} for its \code{tol} argument
}

\examples{
data(ais)
m1 <-  selm(log(Fe) ~ BMI + LBM, family = "sn", data = ais)
intervCP <- confint(m1)
intervDP <- confint(m1, param.type="DP")
confint(m1, parm=2:3)    
confint(m1, parm=c("omega", "alpha"), param.type="DP")                                 
}

