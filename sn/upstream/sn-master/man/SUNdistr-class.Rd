%  file sn/man/SUNdistr-class.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013-2021 Adelchi Azzalini
%---------------------
\name{SUNdistr-class}
\Rdversion{1.1}
\docType{class}
\encoding{UTF-8}

\alias{SUNdistr-class}
\alias{show.SUNdistr}
\alias{show,SUNdistr-method}
\alias{mean.SUNdistr}
\alias{mean,SUNdistr-method}
\alias{vcov.SUNdistr}
\alias{vcov,SUNdistr-method}

\title{Class \code{"SUNdistr"} and its methods}

\description{A class of objects representing Unified Skew-Normal (\acronym{SUN}) 
  distributions.}

\section{Objects from the class}{
  Objects can be created by a call to the function \code{\link{makeSUNdistr}} 
  or by a suitable transformation of some object of this class.}

\section{Slots}{
  \describe{
    \item{\code{dp}:}{a list of parameters of length five, 
      as described in \code{\link{SUNdistr-base}} }
    \item{\code{name}:}{a character string with the name of the multivariate
      variable; it can be an empty string.}
    \item{\code{compNames}:}{a vector of character strings with the names of 
      the component variables.}
   \item{HcompNames}{a vector of character strings with the names of the 
      hidden variables.}    
    }
  }
    
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "SUNdistr-class")}: \dots }
    \item{plot}{\code{signature(x = "SUNdistr-class")}: \dots }
    \item{summary}{\code{signature(object = "SUNdistr-class")}: \dots }
    \item{mean}{\code{signature(x = "SUNdistr")}: \dots}
    \item{vcov}{\code{signature(object = "SUNdistr")}: \dots}
    }
}


\details{See \code{\link{SUNdistr-base}} for a  description of the 
  required structure of \code{dp}.
  
  Note that here the methods \code{mean} and \code{vcov}  are not applied
  to data or to a fitted model, but to a \emph{probability
  distribution},  of which they provide the mean (expected) value
  and the variance-covariance matrix.
  
  The object of this class follow the S4 protocol.
  }

\author{Adelchi Azzalini}

\seealso{
  \code{\link{plot,SUNdistr-method}}, \eqn{\quad}
  \code{\link{summary,SUNdistr-method}}, \eqn{\quad}
  \code{\link{affineTransSUNdistr}},  \code{\link{marginalSUNdistr}}
  
  \code{\link{convertSN2SUNdistr}} to convert a \code{SECdistr} object with
  family \code{"SN"} or \code{"ESN"} to the equivalent \code{SUNdistr-class}
  object
  }

\examples{
  xi <- c(1, 0, -1)
  Omega <- matrix(c(2,1,1, 1,3,1, 1,1,4), 3, 3)
  Delta <- matrix(c(0.72,0.20, 0.51,0.42, 0.88, 0.94), 3, 2, byrow=TRUE)
  Gamma <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
  dp3 <- list(xi=xi, Omega=Omega, Delta=Delta, tau=c(-0.5, 0), Gamma=Gamma)
  sun3 <- makeSUNdistr(dp=dp3, name="firstSUN", compNames=c("x", "w", "z"))
  show(sun3)
  plot(sun3)
  mean(sun3)  # the mean value of the probability distribution
  vcov(sun3)  # the variance-covariance matrix of the probability distribution
  summary(sun3)  # a more detailed summary
}
%----------------
\keyword{classes}
\keyword{distribution}
\keyword{multivariate}
%
\concept{SUN distribution}
\concept{Unified Skew-Normal distribution}
\concept{CSN distribution}
\concept{Closed Skew-Normal distribution}
\concept{FUSN distribution}
\concept{Fundamental Skew-Normal distribution}
