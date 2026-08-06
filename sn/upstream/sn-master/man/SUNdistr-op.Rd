%  file sn/man/SUNdistr-op.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013-2021 Adelchi Azzalini
%---------------------
\name{SUNdistr-op}
\docType{package}
\encoding{UTF-8}

\alias{SUNdistr-op}
\alias{affineTransSUNdistr}
\alias{conditionalSUNdistr}
\alias{convolutionSUNdistr}
\alias{joinSUNdistr}
\alias{marginalSUNdistr}

\title{Operations on SUNdistr-class objects}

\description{
Given an object of \code{SUNdistr-class}, or possibly two such things
in some cases, the functions   perform various operations,
and produce a new object of the same class.}

\usage{
affineTransSUNdistr(object, a, A, name, compNames, HcompNames, drop = TRUE)
conditionalSUNdistr(object, comp, values, eventType = "=", name, drop = TRUE) 
convolutionSUNdistr(object1, object2, name, compNames, HcompNames) 
joinSUNdistr(object1, object2, name, compNames, HcompNames) 
marginalSUNdistr(object, comp, name, drop=TRUE)
}

\arguments{
  \item{object, object1, object2}{objects of class \code{SUNdistr}}
  \item{a}{a numeric vector; see \sQuote{Details}}
  \item{A}{a numeric matrix; see \sQuote{Details}}
  \item{name}{an optional character string with the name of the returned
     distribution}
  \item{compNames}{an optional vector of character strings with the names 
     of the component variables of the returned distribution}
  \item{HcompNames}{an optional vector of character strings with the names 
     of the hidden variables of the returned distribution}
  \item{drop}{a logical value (default: \code{TRUE})  relevant only in the
      case \code{m=1}. When both \code{m=1} and \code{drop=TRUE}, the 
      returned object is of class either \code{SECdistrUv} or \code{SECdistrMv},
      depending on the dimension of the returned object, and family 
      \code{"SN"} or \code{"ESN"}, as appropriate.}
  \item{comp}{a vector of integers representing the selected components}
  \item{values}{a numeric vector which identifies the conditioning event}
  \item{eventType}{a single character value which indicates the type of the
     conditioning event, as described in the \sQuote{Details} section;
     possible values are \code{"="} (default) and \code{">"}}
}

\details{
For  an \code{object} which represents the distribution of a multivariate 
\acronym{SUN} random variable \eqn{Y} of dimension \code{d}, say, a number of
operations are possible, producing a new object of the same class.
This \code{object} could have been created by \code{\link{makeSUNdistr}} 
or it could be the outcome from some previous call to one of the functions
described here.

The function \code{affineTransSUNdistr} computes  the distribution of
\eqn{a+A'Y},  provided \code{A} is a full-rank matrix with
\code{nrow(A)=d} and \code{length(a)=ncol(A)}.
See equation (7.6) of Azzalini & Capitanio (2014).

The function \code{marginalSUNdistr}  builds a \acronym{SUN} distribution 
from the components selected by the \code{comp} vector.

A conditional distribution can be computed using \code{conditionalSUNdistr} 
for two type of events, selected by \code{eventType}. 
The \code{"="} case corresponds to the event \eqn{Y_1=y_1}{Y₁=y₁} where 
\eqn{Y_1}{Y₁} is the subset of components identified
by  the \code{comp} argument, \eqn{y_1}{y₁} is vector specified by the 
\code{values} argument and the equality sign must hold for each component.
See equation (7.6) of Azzalini & Capitanio (2014).

If \code{conditionalSUNdistr} is used with \code{eventType=">"},
the conditiong refers to the event \eqn{Y_1>y_1}{Y₁>y₁}, 
where the inequality must be interpreted components-wise;
see Arellano-Valle & Azzalini (2021) for the underlying mathematical result.
If the conditional distribution is required for the reverse inequality 
condition, \code{"<"} say, 
this is equivalent to consideration of the event \eqn{-Y_1>-y_1}{-Y₁>-y₁}. 
The corresponding distribution can be obtained in two steps: 
first a new variable is constructed reversing the sign of the required
components using \code{affineTransSUNdistr};
then \code{conditionalSUNdistr} is applied to this new variable with
the \code{">"} condition and values \eqn{-y_1}{-y₁}. 
More complex conditions, where the \code{"<"} and \code{">"} signs 
are mixed for different component varables, can be handled similarly, 
by introducing a square matrix \code{A} for \code{affineTransSUNdistr} 
having an appropriate combination of \code{1}s' and \code{-1}'s on its main
diagonal, and 0's elsewhere, and  matching changes of sign to the components  
of \eqn{y_1}{y₁}.

Functions \code{convolutionSUNdistr} and \code{joinSUNdistr} operate under 
the assumptions that \code{object1} and \code{object2} refer to independent 
variables. 
Specifically, \code{convolutionSUNdistr} computes the convolution of the
two objects (i.e. the distribution of the sum of two independent variables),
which must have the same dimension \code{d}. 
Function \code{joinSUNdistr} combines two objects into a joint distribution.

If the arguments \code{name}, \code{compNames} and \code{HcompNames}
are missing, they are composed from the supplied arguments.
}

\value{an object of \code{SUNdistr-class}}

\references{
Arellano-Valle, R. B. and Azzalini, A. (2021).
Some properties of the unified skew-normal distribution.
\emph{Statistical Papers}, 
\doi{https://doi.org/10.1007/s00362-021-01235-2}
and \href{https://arxiv.org/abs/2011.06316}{arXiv:2011.06316}  
  
Azzalini, A. with the collaboration of Capitanio, A. (2014). 
\emph{The Skew-Normal and Related Families}. 
Cambridge University Press, IMS Monographs series.
}

\author{Adelchi Azzalini}

\note{
The present structure and user interface of this function, and of other ones 
related to the \acronym{SUN} distribution, must be considered experimental, 
and they might possibly change in the future.} 

\seealso{
  \code{\link{SUNdistr-base}}, \code{\link{makeSUNdistr}},
  \code{\link{SUNdistr-class}}
}

\examples{
xi <- c(1, 0, -1)
Omega <- matrix(c(2,1,1, 1,3,1, 1,1,4), 3, 3)
Delta <- matrix(c(0.72,0.20, 0.51,0.42, 0.88, 0.94), 3, 2, byrow=TRUE)
Gamma <- matrix(c(1, 0.8, 0.8, 1), 2, 2)
dp3 <- list(xi=xi, Omega=Omega, Delta=Delta, tau=c(-0.5, 0), Gamma=Gamma)
sun3 <- makeSUNdistr(dp=dp3, name="SUN3", compNames=c("x", "w", "z"))
#
a <- c(1,-2)
A <- matrix(1:6, 3, 2)
sun2at <- affineTransSUNdistr(sun3, a, A, "SUN2at", compNames=c("at1", "at2"))
sun2m <- marginalSUNdistr(sun3, comp=c(1,3), name="SUN2m")
sun1c <- conditionalSUNdistr(sun3, comp=c(1,3), values=c(1.1, 0.8),
            eventType=">", name="SUN1c", drop=FALSE)
#
Omega <- matrix(c(5, 1, 1, 6), 2, 2)
Delta <- matrix(c(0.30, 0.50, 0.50, 0.85), 2, 2, byrow=TRUE)
Gamma <- matrix(c(1, 0.18, 0.18, 1), 2, 2)
tau <- c(0.4, -0.8)
dp2 <- list(x=c(1, 0), Omega=Omega, Delta=Delta, tau=tau, Gamma=Gamma)
sun2 <-  makeSUNdistr(dp=dp2, name="SUN2", compNames=c("u", "v"))
#
sun2conv <- convolutionSUNdistr(sun2, sun2m, name="SUN2sum")
sun5 <- joinSUNdistr(sun3, sun2)
}

\keyword{distribution}
\keyword{multivariate} 
\concept{SUN distribution}
\concept{Unified Skew-Normal distribution}
\concept{CSN distribution}
\concept{Closed Skew-Normal distribution}
