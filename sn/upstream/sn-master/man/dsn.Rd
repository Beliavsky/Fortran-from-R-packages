%  file sn/man/dsn.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 1998-2013 Adelchi Azzalini
%---------------------
\name{dsn}
\alias{dsn}
\alias{psn}
\alias{qsn}
\alias{rsn}
\concept{skew-normal distribution}

\title{Skew-Normal Distribution}
\description{Density function, distribution function, quantiles and random
 number generation for the skew-normal (\acronym{SN}) and the extended 
 skew-normal (\acronym{ESN}) distribution.}

\usage{
dsn(x, xi=0, omega=1, alpha=0, tau=0, dp=NULL, log=FALSE)
psn(x, xi=0, omega=1, alpha=0, tau=0, dp=NULL, engine, ...)
qsn(p, xi=0, omega=1, alpha=0, tau=0, dp=NULL, tol=1e-8, solver="NR", ...) 
rsn(n=1, xi=0, omega=1, alpha=0, tau=0,  dp=NULL)
}

\arguments{
 \item{x}{vector of quantiles. Missing values (\code{NA}'s) and \code{Inf}'s
   are allowed.}

 \item{p}{vector of probabilities. Missing values (\code{NA}'s) are allowed}

 \item{xi}{vector of location parameters.}

 \item{omega}{vector of  scale parameters;  must be positive.}

 \item{alpha}{vector of slant parameter(s); \code{+/- Inf} is allowed.
   For \code{psn}, it must be of length 1 if \code{engine="T.Owen"}. 
   For \code{qsn}, it must be of length 1.}

 \item{tau}{a single value representing the `hidden mean' parameter 
   of the \acronym{ESN} distribution; \code{tau=0} (default) corresponds to 
   a \acronym{SN} distribution.}

 \item{dp}{a vector of length 3 (in the \acronym{SN} case) or 
   4 (in the \acronym{ESN} case), whose components represent 
   the individual parameters described above. If \code{dp}
   is specified, the individual parameters cannot be set.}
    
 \item{n}{a positive integer representing the sample size.}

 \item{tol}{a scalar value which regulates the accuracy of the result of 
   \code{qsn}, measured on the probability scale.}

 \item{log}{logical flag used in \code{dsn} (default \code{FALSE}).
   When \code{TRUE}, the logarithm of the density values is returned.}

\item{engine}{a character string  which selects the computing engine;
   this is either \code{"T.Owen"} or \code{"biv.nt.prob"}, the latter from 
   package \code{mnormt}. If \code{tau != 0} or \code{length(alpha)>1}, 
   \code{"biv.nt.prob"} must be used. If this argument is missing, a default 
   selection rule is applied.}

\item{solver}{a character string which selects the numerical method used for 
   solving the quantile equation; possible options are \code{"NR"} (default)
   and \code{"RFB"}, described in the \sQuote{Details} section.}
   
\item{...}{ additional parameters passed to \code{T.Owen}}
}

\value{density (\code{dsn}), probability (\code{psn}), quantile (\code{qsn})
  or random sample (\code{rsn}) from the skew-normal distribution with given
  \code{xi}, \code{omega} and \code{alpha} parameters or from the extended
  skew-normal if \code{tau!=0} }

\section{Details}{
Typical usages are
\preformatted{%
dsn(x, xi=0, omega=1, alpha=0, log=FALSE)
dsn(x, dp=, log=FALSE)
psn(x, xi=0, omega=1, alpha=0,  ...)
psn(x, dp=,  ...)
qsn(p, xi=0, omega=1, alpha=0, tol=1e-8, ...)
qsn(x, dp=, ...)
rsn(n=1, xi=0, omega=1, alpha=0)
rsn(x, dp=)
}
\code{psn} and \code{qsn} make use  of function \code{\link{T.Owen}}
 or \code{\link[mnormt:dmt]{biv.nt.prob}}
 
In \code{qsn}, the choice \code{solver="NR"} selects the Newton-Raphson method 
for solving the quantile equation, while option \code{solver="RFB"}
alternates a step of \emph{regula falsi} with one of bisection. 
The \code{"NR"} method is generally more efficient, but \code{"RFB"} is 
occasionally required in some problematic cases.

In version 1.6-2, the random number generation method for \code{rsn} has
changed;  the so-called transformation method  (also referred to as the
\sQuote{additive representation}) has been adopted for all values of \code{tau}.
Also, the code has been modified so that there is this form of consistency:
provided \code{set.seed()} is reset similarly before calls, code like
\code{rsn(5, dp=1:3)} and \code{rsn(10, dp=1:3)}, for instance, will start with 
the same initial values in the longer sequence as in the shorter sequence.
}

\section{Background}{
The family of skew-normal distributions is an extension of the normal
family, via the introdution of a \code{alpha} parameter which regulates
asymmetry; when \code{alpha=0}, the skew-normal distribution reduces to 
the normal one.  The density function of the \acronym{SN} distribution 
in the \sQuote{normalized}  case having \code{xi=0} and \code{omega=1} is 
\eqn{2\phi(x)\Phi(\alpha x)}, if \eqn{\phi} and \eqn{\Phi} denote the
standard normal density and distribution function.
An early discussion of the skew-normal distribution is given by 
Azzalini (1985); see Section 3.3 for the \acronym{ESN} variant, 
up to a slight difference in the parameterization.

An updated exposition is provided in Chapter 2 of Azzalini and 
Capitanio (2014); the \acronym{ESN} variant is presented Section 2.2. 
See Section 2.3 for an historical account.   
A multivariate version of the distribution is examined in Chapter 5.}


\references{
Azzalini, A. (1985).
  A class of distributions which includes the normal ones.
  \emph{Scand. J. Statist.} \bold{12}, 171-178.

Azzalini, A. with the collaboration of Capitanio, A. (2014). 
 \emph{The Skew-Normal and Related Families}. 
 Cambridge University Press, IMS Monographs series.
}

\seealso{
Functions used by \code{psn}:
   \code{\link{T.Owen}}, \code{\link[mnormt:dmt]{biv.nt.prob}}
 
Related distributions: \code{\link{dmsn}}, \code{\link{dst}}, 
   \code{\link{dmst}}
}
\examples{
pdf <- dsn(seq(-3, 3, by=0.1), alpha=3)
cdf <- psn(seq(-3, 3, by=0.1), alpha=3)
q <- qsn(seq(0.1, 0.9, by=0.1), alpha=-2)
r <- rsn(100, 5, 2, 5)
qsn(1/10^(1:4), 0, 1, 5, 3, solver="RFB")  
}
\keyword{distribution}
