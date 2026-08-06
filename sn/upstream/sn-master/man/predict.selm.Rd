%  file sn/man/profile.selm.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2015 Adelchi Azzalini
%---------------------

\name{predict.selm}
\alias{predict.selm}
\alias{predict,selm-method}   
\concept{confidence interval}
\concept{tolerance interval}
\title{Predict method for selm-class objects}

\description{
Predicted values based on a model object produced by \code{selm} with 
univariate response. 
}
\usage{
  \S4method{predict}{selm}(object, newdata, param.type = "CP",
           interval = "none", level = 0.95, na.action = na.pass, ...) 
}

\arguments{
  \item{object}{an object of class \code{selm} as produced by a call to
    function \code{selm} with univariate response.}
  
  \item{newdata}{an optional data frame in which to look for variables with
    which to predict. If omitted, the fitted values are used.}
    
  \item{param.type}{
    a character string with the required parameterization; it must be one of
    \code{"CP"}, \code{"DP"}, \code{"pseudo-CP"}, or possibly their  
    equivalent lowercase.}
 
  \item{interval}{type of interval calculation among \code{"none", "confidence", 
    "prediction"}; it can be abbreviated.} 
 
  \item{level}{tolerance/confidence level (default value is \code{0.95}).}
  
  \item{na.action}{function determining what should be done with missing 
    values in newdata. The default is to predict \code{NA}.}
    
  \item{...}{not used, only there for compatibility reasons.} 
}

\details{
Predicted values are obtained by evaluating the regression function in the 
dataframe \code{newdata} (which defaults to \code{model.frame(object)}).
Setting \code{interval} other than \code{"none"} produces computation of
confidence or prediction (tolerance) intervals at the specified level.

If \code{newdata} is omitted the predictions are based on the data used for 
the fit.

The action taken in case of missing data is regulated by argument
\code{na.action}, along the lines of function \code{\link[stats]{predict.lm}}.

A detailed description of the methodology underlying \code{predict.selm}
is available in the technical note of Azzalini (2016). 
}

\value{
 a vector of predictions  (if \code{interval="none"}) or a matrix of
 predictions and bounds with column names \code{fit}, \code{lwr}, and 
 \code{upr}, if \code{interval} is set.}
  

\author{Adelchi Azzalini}

\references{
 Azzalini, A. (2016). 
 Derivation of various types of intervals from a \code{selm} object.
 Technical note distributed with the documentation of the \R package
 \code{sn}, in file 
 \href{../doc/selm-intervals.pdf}{\code{selm-intervals.pdf}} within section 
 \sQuote{User guide, package vignettes and other documentation}.
}
 
% \note{}
 
\seealso{
\code{\link{selm}}, \code{\link{summary.selm}},

\code{\link{makeSECdistr}} for the \acronym{CP}/\acronym{DP} parameterizations,

\code{\link[stats]{predict.lm}} for usage of \code{na.action}
}

\examples{
data(barolo)
attach(barolo)  
A75 <- (reseller=="A" & volume==75)
detach(barolo)
m3 <- selm(log(price, 10) ~ age, data=barolo[A75,], family="ST")                
}

