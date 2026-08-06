%  file sn/man/summary.selm.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{summary.selm}
\alias{summary.selm}
\alias{summary.mselm}
\alias{summary,selm-method}    
\alias{summary,mselm-method}
\alias{summary.selm-class}
\alias{summary.mselm-class}
\alias{show,summary.selm-method}
\alias{show,summary.mselm-method}
\title{Summarizing \code{selm} fits}

\description{\code{summary} method for class \code{"selm"} and \code{"mselm"}.}

\usage{
\S4method{summary}{selm}(object, param.type = "CP", cov = FALSE, cor = FALSE)

\S4method{summary}{mselm}(object, param.type = "CP", cov = FALSE, cor = FALSE)
}

 
\arguments{
  \item{object}{an object of class \code{"selm"} or \code{"mselm"} as created
    by a call to function \code{selm}.}
    
  \item{param.type}{a character string which indicates the required type of 
    parameter type; possible values are \code{"CP"} (default), \code{"DP"},
    \code{"pseudo-CP"} and their equivalent lower-case  expressions.}
    
  \item{cov}{a logical value, to indicate if an estimate of the variance and 
    covariance matrix of the estimates is required (default: \code{FALSE}).}
    
  \item{cor}{a logical value, to indicate if an estimate of the correlation 
    matrix of the estimates is required (default: \code{FALSE}).}  
}

\value{An S4 object of class \code{summary.selm} with 12 slots.
    \item{\code{call}:}{the calling statement.}
    \item{\code{family}:}{the parametric family of skew-ellitically 
          contoured distributed (\acronym{SEC}) type.}
    \item{\code{logL}:}{the maximized log-likelihood or penalized 
          log-likelihood value}
    \item{\code{method}:}{estimation method (\code{"MLE"} or \code{"MPLE"})}
    \item{\code{param.type}:}{a characer string with the chosen parameter set.}
    \item{\code{param.table}:}{table of parameters, std.errors and z-values}
    \item{\code{fixed.param}:}{a list of fixed parameter values}
    \item{\code{resid}:}{residual values}
    \item{\code{control}:}{a list with control parameters}
    \item{\code{aux}:}{a list of auxiliary quantities}
    \item{\code{size}:}{a numeric vector with various lengths and dimensions}
    \item{\code{boundary}:}{a logical value which indicates whether the 
          estimates are on the boundary of the parameter space}
}

\note{
There are two reasons why the default choice of \code{param.type} is
\code{CP}.  One is the the easier interpretation of cumulant-based quantities
such as mean value, standard deviation, coefficient of skewness. 

The other reason is more technical and applies only to cases when the 
estimate of the slant parameter \eqn{alpha} of the \acronym{SN} distribution 
is close to the origin: standard asymptotic distribution theory of maximum 
likelihood estimates (MLE's) does not apply in this case and the
corresponding standard errors are not trustworthy. 
The problem is especialy severe at
\eqn{\alpha=0} but to some extent propagates to its vicinity. 
If \eqn{d=1}, adoption of \code{CP} leads to MLE's with regular asymptotic 
distribution across the parameter space, including \eqn{\alpha=0}. 
For \eqn{d>1} and \eqn{\alpha=0,} the problem is still unsolved at the 
present time,  which is the reason why \code{selm} issues a warning
message when the MLE is in the vicinity of \eqn{\alpha=0}; 
see \sQuote{Details} of \code{\link{selm}}.  
For background information, see Sections 3.1.4--6 and 5.2.3 of 
Azzalini and Capitanio (2014)  and references therein.  

This problem does not occur with the the \acronym{SC} and the \acronym{ST}
distribution (unless its tail-weight parameter \code{nu} diverges, that is, 
when we are effectively approaching the \code{SN} case). }

\references{
  Azzalini, A. with the collaboration of Capitanio, A. (2014). 
  \emph{The Skew-Normal and Related Families}. 
  Cambridge University Press, IMS Monographs series.
}

\author{Adelchi Azzalini}

\seealso{
  \code{\link{selm}} function,  
  \code{\linkS4class{selm}} (and \code{mselm}) class, 
  \code{\link{plot.selm}}, \code{\link{dp2cp}},
  \code{\link[=profile.selm]{profile}}, 
  \code{\link[=confint.selm]{confint}}, 
  \code{\link[=predict.selm]{predict}}
}

\examples{ 
data(wines, package="sn")
m5 <- selm(acidity ~ phenols + wine, family="SN", data=wines)
summary(m5)
summary(m5, "dp")
s5 <- summary(m5, "dp", cor=TRUE, cov=TRUE)
dp.cor <- slot(s5, "aux")$param.cor
cov2cor(vcov(m5, "dp")) # the same
#
# m6 <- selm(acidity ~ phenols + wine, family="ST", data=wines) # boundary!?
#
m12 <- selm(cbind(acidity, alcohol) ~ phenols + wine,  family="SN", data=wines)
s12 <- summary(m12)
coef(m12, 'dp')
coef(m12, "dp", vector=FALSE)
#
# see other examples at function selm
}
 
\keyword{regression} 
