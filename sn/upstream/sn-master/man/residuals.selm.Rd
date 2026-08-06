%  file sn/man/summary.selm.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013-2014 Adelchi Azzalini
%---------------------
\name{residuals.selm}
\alias{residuals.selm}
\alias{residuals.mselm}
\alias{residuals,selm-method}             
\alias{residuals,mselm-method}
\alias{fitted.selm}
\alias{fitted.mselm}
\alias{fitted,selm-method}               
\alias{fitted,mselm-method}
\concept{residual values}
\concept{fitted values}
 
\title{Residuals and fitted values from \code{selm} fits}

\description{
  \code{residuals} and \code{fitted} methods for classes \code{"selm"} and 
  \code{"mselm"}.}

\usage{
\S4method{residuals}{selm}(object, param.type = "CP", ...)
\S4method{residuals}{mselm}(object, param.type = "CP", ...)
\S4method{fitted}{selm}(object, param.type = "CP", ...)
\S4method{fitted}{mselm}(object, param.type = "CP", ...)
}

 
\arguments{
  \item{object}{an object of class \code{"selm"} or \code{"mselm"} as created
    by a call to function \code{selm}.}
    
  \item{param.type}{a character string which indicates the required type of 
    parameter type; possible values are \code{"CP"} (default), \code{"DP"},
    \code{"pseudo-CP"} and their equivalent lower-case  expressions.}
    
  \item{...}{not used, included for compatibility with the generic method.}  
}

\value{a numeric vector (for \code{selm-class} objects) or a matrix
     (for \code{mselm-class} objects). }

\note{The possible options of \code{param.type} are described in the
documentation of \code{\link{dp2cp}}; their corresponding outcomes 
differ by an additive constant only. With the \code{"CP"} option (that is,
the \sQuote{centred parametrization}), the residuals are centred around 0, 
at least approximately; this is a reason for setting \code{"CP"} as the 
default option. For more information, see  the \sQuote{Note} in the 
documentation of \code{\link{summary.selm}}. }

\references{
  Azzalini, A. with the collaboration of Capitanio, A. (2014). 
  \emph{The Skew-Normal and Related Families}. 
  Cambridge University Press, IMS Monographs series.
}

\author{Adelchi Azzalini}

\seealso{
  \code{\link{dp2cp}}, \code{\link{summary.selm}}, 
  \code{\link{selm}} function, \code{\linkS4class{selm}-class}
} 


\examples{ 
data(wines, package="sn")
m5 <- selm(acidity ~ phenols + wine, family="SN", data=wines)
residuals(m5)
residuals(m5, "dp")
fitted(m5, "dp")
#
m12 <- selm(cbind(acidity, alcohol) ~ phenols + wine,  family="SN", data=wines)
residuals(m12)
#
# see other examples at function selm
}
 
\keyword{regression} 
