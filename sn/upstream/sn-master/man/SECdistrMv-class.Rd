%  file sn/man/SECdistrMv-class.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{SECdistrMv-class}
\Rdversion{1.1}
\docType{class}
\alias{SECdistrMv-class}
\alias{show,SECdistrMv-method}
\alias{mean,SECdistrMv-method}
\alias{vcov,SECdistrMv-method}

\title{Class \code{"SECdistrMv"}}

\description{A class of objects representing multivariate skew-elliptically 
  contoured (\acronym{SEC}) distributions.}

\section{Objects from the Class}{
  Objects can be created by a call to function \code{\link{makeSECdistr}},
  when its argument \code{dp} is a list, or by a suitable transformation of 
  some object of this class. They can also obtained from an object generated
  by  \code{selm} using the function \code{extractSEDdistr}.}

\section{Slots}{
  \describe{
    \item{\code{family}:}{a character string which identifies the parametric
       family; currently, possible values are: \kbd{"SN"}, \kbd{"ESN"}, 
       \kbd{"ST"}, \kbd{"SC"}.}
    \item{\code{dp}:}{a list of parameters; its length depends on 
       the selected  \code{family}.}
    \item{\code{name}:}{a character string with the name of the multivariate
        variable; it can be an empty string.}
    \item{\code{compNames}:}{a vector of character strings with the names of 
       the component variables.}
    }
  }
    
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "SECdistrMv-class")}: \dots }
    \item{plot}{\code{signature(x = "SECdistrMv-class")}: \dots }
    \item{summary}{\code{signature(object = "SECdistrMv-class")}: \dots }
    \item{mean}{\code{signature(x = "SECdistrMv")}: \dots}
    \item{vcov}{\code{signature(object = "SECdistrMv")}: \dots}
    }
}

\author{Adelchi Azzalini}

\note{See \code{\link{makeSECdistr}} for a detailed description of 
  \code{family} and \code{dp}.
  
  Note that  here methods \code{mean} and \code{vcov}  are not applied
  to data or to a fitted model, but to a \emph{probability
  distribution} instead,  of which they provide the mean (vector) value
  and the variance-covariance matrix.  If methods \code{mean} and \code{vcov} 
  are applied to a distribution for which the mean or the variance do
  not exist, a \code{NULL} value is returned and a warning message
  is issued.}

\seealso{
  \code{\linkS4class{SECdistrUv}},  
  \code{\link{plot,SECdistrMv-method}}, 
  \code{\link{summary,SECdistrMv-method}},
  \code{\link{affineTransSECdistr}},  \code{\link{marginalSECdistr}},
  \code{\link{extractSECdistr}}
}
\examples{
  dp0 <- list(xi=1:2, Omega=diag(3:4), alpha=c(3, -5))
  f10 <- makeSECdistr(dp=dp0, family="SN", name="SN-2D", compNames=c("x", "y"))
  show(f10)
  plot(f10)
  summary(f10)
  mean(f10)  # the mean value of the probability distribution
  vcov(f10)  # the variance-covariance matrix of the probability distribution
}
\keyword{classes}
