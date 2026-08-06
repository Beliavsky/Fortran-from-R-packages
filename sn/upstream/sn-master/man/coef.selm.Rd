%  file sn/man/coef.selm.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2014 Adelchi Azzalini
%---------------------
\name{coef.selm}
\alias{coef.selm}
\alias{coef.mselm}
\alias{coef,selm-method}              % declared in selm-class.Rd
\alias{coef,mselm-method}

 
\title{Coefficients of objects created by \code{selm}}

\description{
  \code{coef} method for classes \code{"selm"} and \code{"mselm"}.}

\usage{
\S4method{coef}{selm}(object, param.type = "CP", ...)
\S4method{coef}{mselm}(object, param.type = "CP", vector=TRUE, ...)
}

 
\arguments{
  \item{object}{an object of class \code{"selm"} or \code{"mselm"} as created
    by a call to function \code{selm}.}
    
  \item{param.type}{a character string which indicates the required type of 
    parameter type; possible values are \code{"CP"} (default), \code{"DP"},
    \code{"pseudo-CP"} and their equivalent lower-case  expressions.}
    
  \item{vector}{a logical value (default is \code{TRUE}) which selects a vector
    or a list format of the retuned value}
    
  \item{...}{not used, included for compatibility with the generic method}  
}

\value{a numeric vector or a list
     (the latter only for \code{mselm-class} objects if \code{vector=FALSE}) }

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
coef(m5)
coef(m5, "dp")
#
m12 <- selm(cbind(acidity, alcohol) ~ phenols + wine,  family="SN", data=wines)
coef(m12)
coef(m12, "DP", vector=FALSE)

}
 
\keyword{regression} 
