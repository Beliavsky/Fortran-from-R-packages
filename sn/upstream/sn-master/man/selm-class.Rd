%  file sn/man/selm-class.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013-2014 Adelchi Azzalini
%---------------------
\name{selm-class}
\Rdversion{1.1}
\docType{class}
\alias{selm-class}               
% \alias{coef,selm-method}           % declared in coef.selm.Rd
\alias{logLik,selm-method}
\alias{plot,selm,ANY-method}
\alias{plot,selm,missing-method}
\alias{show,selm-method}
%\alias{fitted,selm-method}          % declared in residuals.selm.Rd
%\alias{residuals,selm-method}       % declared  in residuals.selm.Rd
\alias{vcov,selm-method}
\alias{weights,selm-method}
%\alias{confint,selm-method}         % declared in confint.selm.Rd
%\alias{predict,selm-method}         % declared in predict.selm.Rd
%\alias{profile,selm-method}         % declared in profile.selm.Rd
%
\alias{mselm-class}
% \alias{coef,mselm-method}          % declared in coef.selm.Rd
\alias{logLik,mselm-method}
\alias{plot,mselm,ANY-method}
\alias{plot,mselm,missing-method}
\alias{show,mselm-method}
%\alias{fitted,mselm-method}         % declared in residuals.selm.Rd
%\alias{residuals,mselm-method}      % declared in residuals.selm.Rd
\alias{vcov,mselm-method}
\alias{weights,mselm-method}

\title{Classes \code{selm} and \code{mselm} of objects created by 
function \code{selm}}

\description{A successful call to function \code{selm} creates an object 
of class either \code{selm} or \code{mselm}, whose structure is described 
in section \sQuote{Slots}. 
A set of methods for object of these classes exist, which are listed 
in the section \sQuote{Methods}. 
Specific documentation is provided separately for the methods involving 
a more extended description; see the section \sQuote{See also}. }

\section{Objects from the class}{
An object can be created by a successful call to function \code{selm}.}

\section{Slots}{
   \describe{
    \item{\code{call}:}{the calling statement.}
    \item{\code{family}:}{the parametric family of skew-ellitically 
          contoured distributed (SEC) type.}
    \item{\code{logL}:}{log-likelihood or penalized log-likelihood value
         achieved at the end of the maximization process.}
    \item{\code{method}:}{estimation method (\code{"MLE"} or \code{"MPLE"}).}
    \item{\code{param}:}{estimated parameters, for various parameterizations.}
    \item{\code{param.var}:}{approximate variance matrices of the parameter 
         estimates, for various parameterizations.}
    \item{\code{size}:}{a numeric vector with size of various components.}
    \item{\code{fixed.param}:}{a vector of parameters which have been kept
        fixed in the fitting process, if any.}
    \item{\code{residuals.dp}:}{residual values, for DP-type parameters.}
    \item{\code{fitted.values.dp}:}{fitted values, for DP-type parameters.}
    \item{\code{control}:}{a list with control parameters.}
    \item{\code{input}:}{a list of selected input values.}
    \item{\code{opt.method}:}{a list with details on the numerical optimization 
          method used in the fitting process.}
  }
}
\section{Methods}{
  \tabular{ll}{
    \code{coef} \tab \code{signature(object = "selm")}: ... \cr
    \code{logLik} \tab \code{signature(object = "selm")}: ... \cr
    % \code{plot} \tab \code{signature(x = "selm", y = "ANY")}: ... \cr
    % \code{plot} \tab \code{signature(x = "selm", y = "missing")}: ... \cr
    \code{plot} \tab \code{signature(x = "selm")}: ... \cr
    \code{show} \tab \code{signature(object = "selm")}: ... \cr
    \code{summary} \tab \code{signature(object = "selm")}: ... \cr
    \code{residuals} \tab \code{signature(object = "selm")}: ... \cr
    \code{fitted} \tab \code{signature(object = "selm")}: ... \cr
    \code{vcov} \tab \code{signature(object = "selm")}: ... \cr
    \code{weights} \tab \code{signature(object = "selm")}: ... \cr
    \code{confint} \tab \code{signature(object = "selm")}: ... \cr
    \code{predict} \tab \code{signature(object = "selm")}: ... \cr
    \code{profile} \tab \code{signature(fitted = "selm")}: ... \cr
    \tab \cr
    \code{coef} \tab \code{signature(object = "mselm")}: ... \cr
    \code{logLik} \tab \code{signature(object = "mselm")}: ... \cr
    \code{plot} \tab \code{signature(x = "mselm")}: ... \cr
    \code{show} \tab \code{signature(object = "mselm")}: ... \cr
    \code{summary} \tab \code{signature(object = "mselm")}: ... \cr
    \code{residuals} \tab \code{signature(object = "mselm")}: ... \cr
    \code{fitted} \tab \code{signature(object = "mselm")}: ... \cr
    \code{vcov} \tab \code{signature(object = "mselm")}: ... \cr
    \code{weights} \tab \code{signature(object = "mselm")}: ... \cr
	}
}

%\references{%%  ~~put references to the literature/web site here~~}

\author{Adelchi Azzalini}

\note{See \code{\link{dp2cp}} for a description of possible parameter sets.
When \code{logLik} is used on an object obtained using the MPLE estimation
method, the value reported is actually the \emph{penalized} log-likelihood.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{See also  the  \code{\link{selm}} function,  
  and related methods \code{\link[=plot.selm]{plot}},
  \code{\link[=summary.selm]{summary}}, % \code{\linkS4class{summary.selm}},  
  % \code{\link{dp2cp}},
  \code{\link[=profile.selm]{profile}}, 
  \code{\link[=confint.selm]{confint}}, 
  \code{\link[=predict.selm]{predict}}
}

\examples{
data(ais)
m1 <- selm(log(Fe) ~ BMI + LBM, family="SN", data=ais)
summary(m1)
plot(m1)
logLik(m1)
res <- residuals(m1)
fv <- fitted(m1)
# 
data(wines, package="sn")
m2 <- selm(alcohol ~ malic + phenols, data=wines)
#
m12 <- selm(cbind(acidity, alcohol) ~ phenols + wine, family="SN", data=wines)
coef(m12)
cp <- coef(m12, vector=FALSE)
dp <- coef(m12, "DP", vector=FALSE)
plot(m12)
plot(m12, which=2, col="gray60", pch=20)
}

\keyword{classes}
