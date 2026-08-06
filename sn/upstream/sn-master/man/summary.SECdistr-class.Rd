%  file sn/man/summary.SECdistr-class.Rd
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{summary.SECdistrMv-class}
\Rdversion{1.1}
\docType{class}
\alias{summary.SECdistrMv-class}
\alias{summary.SECdistrUv-class}
\alias{show,summary.SECdistrMv-method}
\alias{show,summary.SECdistrUv-method}

\title{Classes \code{summary.SECdistrMv} and \code{summary.SECdistrUv}}

\description{Summaries of objects of classes \code{SECdistrMv} and 
\code{SECdistrUv}}

\section{Objects from the Class}{
Objects can be created by calls of type \code{summary(object)} when 
\code{object} is  of class either \code{"SECdistrMv"} or \code{"SECdistrUv"}.}

\section{Slots}{
  \describe{
    \item{\code{family}:}{A character string which represents
          the parametric family of \acronym{SEC} type }
    \item{\code{dp}:}{Object of class \code{"list"} or \code{"vector"} for
        \code{"SECdistrMv"} and \code{"SECdistrUv"}, respectively}
    \item{\code{name}:}{Object of class \code{"character"} with the name of 
        distribution }
    \item{\code{compNames}:}{For \code{"SECdistrMv"} objects, a character 
        vector with names of the components of the multivariate distribution}
    \item{\code{cp}:}{Object of class \code{"list"} or \code{"vector"} for
        \code{"SECdistrMv"} and \code{"SECdistrUv"}, respectively}
    \item{\code{cp.type}:}{a character string of the \acronym{CP} version}
    \item{\code{aux}:}{A list of auxiliary quantities }
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "summary.SECdistrMv")}: ... }
    \item{show}{\code{signature(object = "summary.SECdistrUv")}: ... }
	 }
}

%\references{%%  ~~put references to the literature/web site here~~}

\author{Adelchi Azzalini}

% \note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{summary.SECdistrMv}}, \code{\link{summary.SECdistrUv}},

\code{\link{makeSECdistr}}, \code{\link{dp2cp}} }

% \examples{showClass("summary.SECdistrMv")}
\keyword{classes}
