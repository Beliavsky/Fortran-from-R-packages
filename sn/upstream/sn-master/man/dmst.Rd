%  file sn/man/dmst.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2002-2013 Adelchi Azzalini
%---------------------
\name{dmst}
\alias{dmst}
\alias{pmst}
\alias{rmst}
\alias{dmsc}
\alias{pmsc}
\alias{rmsc}
\concept{skew-t distribution}
\concept{multivariate skew-t distribution}

\title{Multivariate skew-\eqn{t} distribution and skew-Cauchy distribution}

\description{Probability density function, distribution function and random 
number generation for the multivariate skew-\eqn{t} (\acronym{ST}) and
skew-Cauchy (\acronym{SC}) distributions.}

\usage{
dmst(x, xi=rep(0,length(alpha)), Omega, alpha, nu=Inf, dp=NULL, log=FALSE)
pmst(x, xi=rep(0,length(alpha)), Omega, alpha, nu=Inf, dp=NULL, ...)
rmst(n=1, xi=rep(0,length(alpha)), Omega, alpha, nu=Inf, dp=NULL)
dmsc(x, xi=rep(0,length(alpha)), Omega, alpha, dp=NULL, log=FALSE)
pmsc(x, xi=rep(0,length(alpha)), Omega, alpha, dp=NULL, ...)
rmsc(n=1, xi=rep(0,length(alpha)), Omega, alpha, dp=NULL)
}

\arguments{

 \item{x}{for \code{dmst} and \code{dmsc}, this is either a vector of length
   \code{d}, where \code{d=length(alpha)}, or a matrix with \code{d} columns,
   representing the coordinates of the point(s) where the density must be
   avaluated; for \code{pmst} and \code{pmsc}, only a vector of length
   \code{d} is allowed.}
   
 \item{xi}{a numeric vector of length \code{d} representing the 
   location parameter of the distribution; see \sQuote{Background}.
   In a call to \code{dmst} or \code{dmsc}, \code{xi} can be a matrix,
   whose rows represent a set of location parameters; 
   in this case, its dimensions must match those of \code{x}.}
  

 \item{Omega}{a symmetric positive-definite matrix of dimension \code{(d,d)}; 
   see Section \sQuote{Background}.}

 \item{alpha}{a numeric vector of length \code{d} which regulates the slant
   of the density; see Section \sQuote{Background}.
   \code{Inf} values in \code{alpha} are not allowed.}

 \item{nu}{a positive value representing the degrees of freedom of
   \acronym{ST} distribution; does not need to be integer. 
   Default value is \code{nu=Inf} which corresponds to the multivariate
   skew-normal distribution.}

 \item{dp}{a list with three elements named \code{xi}, \code{Omega},
  \code{alpha} and \code{nu}, containing quantities as described above. If
  \code{dp} is specified, this prevents specification of the individual
  parameters.}

 \item{n}{a numeric value which represents the number of random vectors to be
  drawn; default value is \code{1}.}

 \item{log}{logical (default value: \code{FALSE}); if \code{TRUE},
  log-densities are returned.}

 \item{...}{additional parameters passed to \code{pmt}; 
     see \sQuote{Details} for further information.}

}

\value{A vector of density values (\code{dmst} and \code{dmsc}) or a single
   probability (\code{pmst} and \code{pmsc}) or a matrix of random points
   (\code{rmst} and \code{rmsc}).}

\details{Typical usages are
\preformatted{%
dmst(x, xi=rep(0,length(alpha)), Omega, alpha, nu=Inf, log=FALSE)
dmst(x, dp=, log=FALSE)
pmst(x, xi=rep(0,length(alpha)), Omega, alpha, nu=Inf, ...)
pmst(x, dp=, ...)
rmst(n=1, xi=rep(0,length(alpha)), Omega, alpha, nu=Inf)
rmst(n=1, dp=)
dmsc(x, xi=rep(0,length(alpha)), Omega, alpha, log=FALSE)
dmsc(x, dp=, log=FALSE)
pmsc(x, xi=rep(0,length(alpha)), Omega, alpha, ...)
pmsc(x, dp=, ...)
rmsc(n=1, xi=rep(0,length(alpha)), Omega, alpha)
rmsc(n=1, dp=)
}
For efficiency reasons, \code{rmst, rmsc} make very limited checks on the 
validity of the arguments. For instance, failure to positive definiteness 
of \code{Omega} would not be detected, and an uncontrolled crash occurs.
Function \code{pmst} requires \code{\link[mnormt]{dmt}} from package 
\pkg{mnormt}; the accuracy of its computation can be controlled via 
argument \code{\dots}.}

\section{Background}{
The family of multivariate \acronym{ST} distributions is an extension of the 
multivariate Student's \eqn{t} family, via the introduction of a \code{alpha} 
parameter which regulates asymmetry; when \code{alpha=0}, the skew-\eqn{t}
distribution reduces to the commonly used form of multivariate Student's 
\eqn{t}.  Further,  location is regulated by \code{xi} and scale by 
\code{Omega}, when its diagonal terms are not all 1's.
When \code{nu=Inf} the distribution reduces to the multivariate skew-normal 
one; see \code{dmsn}. Notice that the location vector \code{xi}
does not represent the mean vector of the distribution (which in fact
may not even exist if \code{nu <= 1}), and similarly \code{Omega} is not 
\emph{the} covariance matrix of the distribution, although it is \emph{a} 
covariance matrix. 
For additional information, see Section 6.2 of the monograph indicated 
in the reference list.

The family of multivariate \acronym{SC} distributions is the subset of the
\acronym{ST} family, obtained when \code{nu=1}. While in the univariate case
there are specialized functions for the \acronym{SC} distribution,
\code{dmsc}, \code{pmsc} and \code{rmsc} simply make a call to \code{dmst,
pmst, rmst} with argument \code{nu} set equal to 1.}

\section{Details}{
Computation of  \code{pmst} can be regulated via \code{...}, by tuning the
arguments which control numerical accuracy  \code{pmt} (for the case of 
integer \code{nu}) or those of \code{integrate} 
(for the case of fractional \code{nu}).
Computation of \code{pmsc} is performed by \code{pmst} setting \code{nu=1}.
}
\references{
%  Azzalini, A. and Capitanio, A. (2003).
%  Distributions generated by perturbation of symmetry 
%  with emphasis on a multivariate skew \emph{t} distribution.
%  \emph{J.Roy. Statist. Soc. B} \bold{65}, 367--389.
  
  Azzalini, A. with the collaboration of Capitanio, A. (2014). 
  \emph{The Skew-Normal and Related Families}. 
  Cambridge University Press, IMS Monograph series.
}
\seealso{
\code{\link{dst}},  \code{\link{dsc}}, \code{\link{dmsn}}, 
\code{\link[mnormt]{dmt}}, \code{\link[mnormt]{pmt}}, 
\code{\link{makeSECdistr}}, \code{\link[stats]{integrate}}
}
\examples{
x <- seq(-4,4,length=15)
xi <- c(0.5, -1)
Omega <- diag(2)
Omega[2,1] <- Omega[1,2] <- 0.5
alpha <- c(2,2)
pdf <- dmst(cbind(x,2*x-1), xi, Omega, alpha, 5)
rnd <- rmst(10,  xi, Omega, alpha, 6)
p1 <- pmst(c(2,1), xi, Omega, alpha, nu=5)
p2 <- pmst(c(2,1), xi, Omega, alpha, nu=5, abseps=1e-12, maxpts=10000)
}
\keyword{distribution}
\keyword{multivariate}

