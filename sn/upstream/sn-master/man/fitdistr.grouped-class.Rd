%  file sn/man/fitdistr.grouped-class.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2022 Adelchi Azzalini
%---------------------
\name{fitdistr.grouped-class}
\docType{class}
\alias{fitdistr.grouped-class}  
\alias{logLik.fitdistr.grouped}  
%\alias{logLik,fitdistr.grouped-method}
\alias{coef.fitdistr.grouped}         
% \alias{coef,fitdistr.grouped-method} 
\alias{vcov.fitdistr.grouped}   
%\alias{vcov,fitdistr.grouped-method}        
%%\alias{plot.fitdistr.grouped}  % documented separately
\alias{print.fitdistr.grouped}  
%\alias{print,fitdistr.grouped-method}
\alias{summary.fitdistr.grouped}  
\alias{fitted.fitdistr.grouped}  
% \alias{fitted,fitdistr.grouped-method}      

%
\title{Methods for objects of class created by \code{fitdistr.grouped}}

\description{A successful call to function \code{fitdistr.grouped} creates 
an object of class, also named \code{fitdistr.grouped}, for which a set 
of methods exist.
The structure of an object of this class is described in section 
\sQuote{Object components}.
}

\usage{
  \S3method{logLik}{fitdistr.grouped}(object, ...)
  \S3method{coef}{fitdistr.grouped}(object, ...)
  \S3method{vcov}{fitdistr.grouped}(object, ...)
  \S3method{print}{fitdistr.grouped}(x, ...)
  \S3method{summary}{fitdistr.grouped}(object, cor=FALSE, ...)
  \S3method{fitted}{fitdistr.grouped}(object, full=FALSE, ...)
}

\arguments{
  \item{x, object}{an object of class \code{fitdistr.grouped} as created 
     by a call to the function with this name.}
  \item{cor}{logical (default=\code{FALSE}); is the correlation matrix required?}   
  \item{full}{logical (default=\code{FALSE}); must the vector of fitted 
       frequencies include the boundary classes, when they are added to 
       cover the full support of the fitted distribution?}
   \item{...}{further arguments passed to or from other methods.}     
  }

\section{Object components}{
  Components of an object of class \code{fitdistr.grouped}:
  \describe{
  \item{\code{call}}{the matched call}
  \item{\code{family}}{the selected \code{family} of distributions}
  \item{\code{logL}}{the achieved maximum log-likelihood}
  \item{\code{param}}{a vector of estimated parameters}
  \item{\code{vcov}}{the approximate variance-covariance matrix of the estimates}
  \item{\code{input}}{a list with the input quantities and some derived ones}
  \item{\code{opt}}{a list as returned by \code{optim}}
  }
}
\author{Adelchi Azzalini}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 the function \code{\link{fitdistr.grouped}}, 
 the plotting method \code{\link{plot.fitdistr.grouped}}
 }

\references{
Possolo, A., Merkatas, C. and Bodnar, O. (2019).
Asymmetrical uncertainties.
\emph{Metrologia} 56, 045009.
}

\examples{
data(barolo)
attach(barolo)
A75 <- (reseller=="A" & volume==75)
logPrice <- log(price[A75], 10) # as used in selm documentation; see its fitting
detach(barolo)
breaks <- seq(1, 3, by=0.25)
f <- cut(logPrice, breaks = breaks)
counts <- tabulate(f, length(levels(f))) 
fit.logPrice.gr <- fitdistr.grouped(breaks, counts, family='ST')
summary(fit.logPrice.gr) # compare this fit with the ungrouped data fitting 
print(fit.logPrice.gr)
coef(fit.logPrice.gr)
vcov(fit.logPrice.gr)
data.frame(intervals=levels(f), counts, fitted=format(fitted(fit.logPrice.gr)))
full.intervals <- c("(-Inf, 1]", levels(f), "(3, Inf)")
data.frame("full-range intervals" = full.intervals,
       "full-range counts" = c(0, counts, 0), 
      "full-range fit" = fitted(fit.logPrice.gr, full=TRUE))
sum(counts) - sum(fitted(fit.logPrice.gr))  
sum(counts) - sum(fitted(fit.logPrice.gr, full=TRUE))  # must be "nearly 0" 
#---
# Use first entry in Table 3 of Possolo et al. (2019) and do similar fitting
# to the *probability* values, not observation counts
afcrc59 <- 1.141
breaks <- c(-Inf, afcrc59 - 0.033, afcrc59, afcrc59 + 0.037, Inf)
prob <-  c(0.16, 0.50, 0.84) 
cum.percent <- c(0, prob, 1)*100 
fitSN <- fitdistr.grouped(breaks, counts=diff(cum.percent), family="SN") 
print(coef(fitSN))
print(rbind(target=c(prob, 1)*100, fitted=cumsum(fitted(fitSN))), digits=5)
# Note: given the nature of these data (i.e. probabilities, not counts), there
#       is no point in using vcov, logLik and summary on the fitted object.
}

\keyword{classes}
\concept{grouped data}
