%  file sn/man/summary.SUNdistr-class.Rd
%  This file is a component of the package 'sn' for R
%  copyright (C) 2021 Adelchi Azzalini
%---------------------
\name{summary.SUNdistr-class}

\Rdversion{1.1}
\docType{class}
\encoding{UTF-8}

\alias{summary.SUNdistr-class}
\alias{show,summary.SUNdistr-method}

\title{Class \code{summary.SUNdistr}}

\description{Summaries of objects of classes \code{SUNdistr}}

\section{Objects from the Class}{
Objects can be created by calls of type \code{summary(object)} when 
\code{object} is  of class \code{"SUNdistr"}.}

\section{Slots}{
 \describe{
  \item{\code{dp}:}{a list of parameters}
  \item{name}{the name of the distribution, a character string}
  \item{compNames}{the names of the components, a character vector}
  \item{HcompNames}{the names of the hidden components, a character vector}   
  \item{mean}{the mean value, a vector}
  \item{var.cov}{the variance-covariance matrix}
  \item{gamma1}{the marginal indices of asymmetry, a vector}
  \item{cum3}{the third order cumulants, a three-dimensional array}
  \item{mardia}{the Mardia's measures of multivariate asymmetry and
     skewness, a vector of length two}
  }
}

\section{Methods}{
  \describe{ 
     \item{show}{\code{signature(object = "summary.SUNdistr")}: ... }
	 }
}

%\references{%%  ~~put references to the literature/web site here~~}

\author{Adelchi Azzalini}

% \note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link{summary.SUNdistr}}, \code{\link{makeSUNdistr}}}

% \examples{showClass("summary.SUNdistrMv")}
\keyword{classes}
