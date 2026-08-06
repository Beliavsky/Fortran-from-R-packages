%  file sn/man/dmsn.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 1998-2013 Adelchi Azzalini
%---------------------
\name{dmsn}
\alias{dmsn}
\alias{pmsn}
\alias{rmsn}
\concept{skew-normal distribution}
\concept{multicariate skew-normal distribution}

\title{Multivariate skew-normal distribution}
\description{
Probability density function, distribution function and random number 
generation for the multivariate skew-normal (\acronym{SN}) distribution.
}

\usage{
dmsn(x, xi=rep(0,length(alpha)), Omega, alpha, tau=0, dp=NULL, log=FALSE)
pmsn(x, xi=rep(0,length(alpha)), Omega, alpha,  tau=0, dp=NULL, ...)
rmsn(n=1, xi=rep(0,length(alpha)), Omega, alpha,  tau=0, dp=NULL)
}

\arguments{

 \item{x}{either a vector of length \code{d}, 
   where \code{d=length(alpha)}, or a matrix with \code{d} columns, 
   giving the coordinates of the point(s) where the density or the
   distribution function must be evaluated.}

 \item{xi}{a numeric vector of length \code{d} representing the 
   location parameter of the distribution; see \sQuote{Background}.
   In a call to \code{dmsn} and \code{pmsn}, \code{xi} can be a matrix,
   whose rows represent a set of location parameters;
   in this case, its dimensions must match those of \code{x}.}

 \item{Omega}{a symmetric positive-definite matrix of dimension \code{(d,d)};
   see \sQuote{Background}.}

 \item{alpha}{a numeric vector which regulates the slant of the density; 
   see \sQuote{Background}. \code{Inf} values in \code{alpha} are not allowed.}

 \item{tau}{a single value representing the `hidden mean' parameter 
   of the \acronym{ESN} distribution; \code{tau=0} (default) corresponds to 
   a \acronym{SN} distribution.}

 \item{dp}{
   a list with three elements, corresponding to \code{xi}, \code{Omega} and 
   \code{alpha} described above; default value \code{FALSE}. 
   If \code{dp} is assigned, individual parameters must not be specified. }

 \item{n}{a numeric value which represents the number of random vectors
   to be drawn.}

\item{log}{logical (default value: \code{FALSE}); 
  if \code{TRUE}, log-densities are returned.}

\item{...}{additional parameters passed to \code{\link[mnormt]{pmnorm}}.}

} % end arguments

\value{
A vector of density values (\code{dmsn}) or of probabilities
(\code{pmsn}) or a matrix of random  points (\code{rmsn}).
}
\details{Typical usages are
\preformatted{%
dmsn(x, xi=rep(0,length(alpha)), Omega, alpha, log=FALSE)
dmsn(x, dp=, log=FALSE)
pmsn(x, xi=rep(0,length(alpha)), Omega, alpha, ...)
pmsn(x, dp=)
rmsn(n=1, xi=rep(0,length(alpha)), Omega, alpha)
rmsn(n=1, dp=)
}
For efficiency reasons, \code{rmsn} makes very limited checks on the 
validity of the arguments. For instance, failure to positive definiteness 
of \code{Omega} would not be detected, and an uncontrolled crash occurs. 
Function \code{pmsn} makes use of \code{pmnorm} from package \pkg{mnormt};
the accuracy of its computation can be controlled via \code{...}
}

\section{Background}{
The multivariate skew-normal distribution is discussed by Azzalini and 
Dalla Valle (1996). The \code{(Omega,alpha)}
parametrization adopted here is the one of Azzalini and Capitanio (1999).
Chapter 5 of Azzalini and Capitanio (2014) provides an extensive account,
including subsequent developments.

Notice that the location vector \code{xi} does not represent the mean vector 
of the distribution. Similarly, \code{Omega} is not \emph{the} covariance 
matrix of the distribution, although it is \emph{a} covariance matrix. 
Finally, the components of \code{alpha} are not equal to the slant parameters 
of the marginal  distributions; to fix the marginal parameters at prescribed 
values, it is convenient to start from the OP parameterization, as illustrated 
in the \sQuote{Examples} below. Another option is to start from the \acronym{CP}
parameterization, but notice that, at variance from the \acronym{OP}, not all
\acronym{CP} sets are invertible to lend a \acronym{DP} set.}

\references{

Azzalini, A. and Capitanio, A. (1999).
Statistical applications of the multivariate skew normal distribution.
\emph{J.Roy.Statist.Soc. B} \bold{61}, 579--602. Full-length version 
available at \url{https://arXiv.org/abs/0911.2093}

Azzalini, A. with the collaboration of Capitanio, A. (2014). 
\emph{The Skew-Normal and Related Families}. 
Cambridge University Press, IMS Monographs series.

Azzalini, A. and Dalla Valle, A. (1996).
The multivariate skew-normal distribution.
\emph{Biometrika} \bold{83}, 715--726.

}

\seealso{\code{\link{dsn}},  \code{\link{dmst}},  \code{\link[mnormt]{pmnorm}},
         \code{\link{op2dp}},  \code{\link{cp2dp}}}

\examples{
x <- seq(-3,3,length=15)
xi <- c(0.5, -1)
Omega <- diag(2)
Omega[2,1] <- Omega[1,2] <- 0.5
alpha <- c(2,-6)
pdf <- dmsn(cbind(x, 2*x-1), xi, Omega, alpha)
cdf <- pmsn(cbind(x, 2*x-1), xi, Omega, alpha)
p1 <- pmsn(c(2,1), xi, Omega, alpha)
p2 <- pmsn(c(2,1), xi, Omega, alpha, abseps=1e-12, maxpts=10000)
#
rnd <- rmsn(10, xi, Omega, alpha)
#
# use OP parameters to fix marginal shapes at given lambda values:
op <- list(xi=c(0,1), Psi=matrix(c(2,2,2,3), 2, 2), lambda=c(5, -2))
rnd <- rmsn(10, dp=op2dp(op,"SN"))
# 
# use CP parameters to fix mean vector, variance matrix and marginal skewness:
cp <- list(mean=c(0,0), var.cov=matrix(c(3,2,2,3)/3, 2, 2), gamma1=c(0.8, 0.4))
dp <- cp2dp(cp, "SN")
rnd <- rmsn(5, dp=dp)
}
\keyword{distribution}
\keyword{multivariate}
