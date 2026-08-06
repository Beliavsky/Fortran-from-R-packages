%  file sn/man/makeSECdistr.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013-2014 Adelchi Azzalini
%---------------------
\name{makeSECdistr}
\encoding{UTF-8}
\alias{makeSECdistr}
\concept{skew-elliptical distribution}
\title{Build a skew-elliptically contoured distribution}

\description{Build an object which identifies a skew-elliptically contoured 
  distribution (\acronym{SEC}), in the univariate and in the multivariate case.
  The term \sQuote{skew-elliptical distribution} is a synonym of \acronym{SEC} 
  distribution.}

\usage{makeSECdistr(dp, family, name, compNames)}

\arguments{

  \item{dp}{a numeric vector (in the univariate case) or a list (in the
    multivariate case) of parameters which identify the specific distribution
    within the named \code{family}. See \sQuote{Details} for their expected 
    structure.}

  \item{family}{a character string which identifies the parametric
    family; currently, possible values are: \kbd{"SN"}, \kbd{"ESN"}, 
    \kbd{"ST"}, \kbd{"SC"}. 
    See \sQuote{Details} for additional information.}

  \item{name}{an optional character string with the name of the distribution.
    If missing, one is created.}

  \item{compNames}{in the multivariate case, an optional vector of character
    strings with the names of the component variables; its length must be 
    equal to the dimensionality of the distribution being generated. 
    If missing and the first component of \code{dp} is a named vector,
    its names are used as \code{compNames}; otherwise
    the components are named \code{"V1"}, \code{"V2"}, \dots}
}

\details{If \code{dp} is a numeric vector, a univariate distribution is built. 
  Alternatively, if \code{dp} is a list, a multivariate distribution is 
  built. In both cases, the required number of components of \code{dp} 
  depends on \code{family}: it must be  \code{3} for \kbd{"SN"} and 
  \kbd{"SC"}; it must be \code{4} for \kbd{"ESN"} and \kbd{"ST"}.

  In the univariate case, the first three components of \code{dp} represent 
  what for the specific distributions are denoted \code{xi} (location), 
  \code{omega} (scale, positive) and \code{alpha} (slant); see functions 
  \code{\link{dsn}}, \code{\link{dst}}, \code{\link{dsc}} for their
  description. 
  The fourth component, when it exists, represents either \code{tau} 
  (hidden variable mean) for \kbd{"ESN"} or \code{nu} (degrees of freedom) 
  for \kbd{"ST"}. The names of the individual parameters are attached
  to the components of \code{dp} in the returned object.
  
  In the multivariate case, \code{dp} is a list with components having 
  similar role as in the univariate case, but \code{xi=dp[[1]]} and 
  \code{alpha=dp[[3]]} are now vectors and the scale parameter
  \code{Omega=dp[[2]]} is a symmetric positive-definite matrix. 
  For a multivariate distribution of dimension 1 (which can be created, 
  although a warning message is issued), \code{Omega} corresponds to the 
  square of \code{omega} in the univariate case. 
  Vectors \code{xi} and \code{alpha} must be of length \code{ncol(Omega)}.
  See also functions \code{\link{dmsn}}, \code{\link{dmst}} and 
  \code{\link{dmsc}}.
  The fourth component, when it exists, is a scalar with the same role as 
  in the univariate case. 
  
  In the univariate case \code{alpha=Inf} is allowed, but in the multivariate
  case all components of the vector \code{alpha} must be finite.
  
  An object built by this function operates according to the S4 protocol.
}

\section{Background}{
  For background information, see Azzalini and Capitanio (2014), specifically
  Chapters 2 and 4 for univariate cases, Chapters 5 and 6 for multivariate
  cases; Section 6.1 provides a general formulation of \acronym{SEC}
  distributions.
  
  If the slant parameter \code{alpha} is \code{0} (or a vector of \code{0}'s,
  in the multivariate case), the distribution is of classical elliptical 
  type.
  
  The \acronym{ESN} distribution is included here  as a members of the
  \acronym{SEC} class, with a very slight extension of the original definition
  of this class, since the only difference is the non-zero truncation point
  of the unobserved component of the \code{(d+1)}-dimensional \acronym{EC} 
  variable.
  }
  
\value{In the univariate case, an object of class \code{SECdistrUv};
  in the multivariate case, an object of class \code{SECdistrMv}.
  See \code{\link{SECdistrUv-class}} and \code{\link{SECdistrMv-class}}
  for their description.
}

\references{
  Azzalini, A. with the collaboration of Capitanio, A. (2014). 
  \emph{The Skew-Normal and Related Families}. 
  Cambridge University Press, IMS Monographs series.
}

\author{Adelchi Azzalini}

\seealso{

 The description of classes  \code{\link{SECdistrUv-class}} and
 \code{\link{SECdistrMv-class}}

 \code{\link{plot.SECdistr}} for plotting and 
 \code{\link{summary.SECdistr}} for summaries

 Related functions \code{\link{dsn}}, \code{\link{dst}}, \code{\link{dsc}}, 
 \code{\link{dmsn}}, \code{\link{dmst}},  \code{\link{dp2cp}}
 
 Functions \code{\link{affineTransSECdistr}} and
  \code{\link{conditionalSECdistr}} to manipulate objects of class    
  \code{\link{SECdistrMv-class}} 
 
 Function \code{\link{extractSECdistr}} to extract objects of class 
 \code{\link{SECdistrMv-class}} and \code{\link{SECdistrUv-class}} 
 representing the \acronym{SEC} distribution of a \code{\link{selm}} fit
 
}

\examples{    
f1 <- makeSECdistr(dp=c(3,2,5), family="SN", name="First-SN")
show(f1)
summary(f1)
plot(f1)
plot(f1, probs=c(0.1, 0.9))
#
f2 <- makeSECdistr(dp=c(3, 5, -4, 8), family="ST", name="First-ST")
f9 <- makeSECdistr(dp=c(5, 1, Inf, 0.5), family="ESN", name="ESN,alpha=Inf")
#
dp0 <- list(xi=1:2, Omega=diag(3:4), alpha=c(3, -5))
f10 <- makeSECdistr(dp=dp0, family="SN", name="SN-2d", compNames=c("u1", "u2"))
#
dp1 <- list(xi=1:2, Omega=diag(1:2)+outer(c(3,3),c(2,2)), alpha=c(-3, 5), nu=6)
f11 <- makeSECdistr(dp=dp1, family="ST", name="ST-2d", compNames=c("t1", "t2"))
} 

\keyword{distribution}
\keyword{multivariate} 
