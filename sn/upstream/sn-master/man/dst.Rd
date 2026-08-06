%  file sn/man/dst.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2002-2023 Adelchi Azzalini
%---------------------
\name{dst}
\alias{dst}
\alias{pst}
\alias{qst}
\alias{rst}
\concept{skew-t distribution}

\title{Skew-\eqn{t} Distribution}
\description{Density function, distribution function, quantiles and 
  random number generation for the skew-\eqn{t} (\acronym{ST}) distribution.}

\usage{
dst(x, xi=0, omega=1, alpha=0, nu=Inf, dp=NULL, log=FALSE) 
pst(x, xi=0, omega=1, alpha=0, nu=Inf, dp=NULL, method=0, lower.tail=TRUE, 
    log.p=FALSE, ...)
qst(p, xi=0, omega=1, alpha=0, nu=Inf, tol=1e-08, dp=NULL, method=0, ...)
rst(n=1, xi=0, omega=1, alpha=0, nu=Inf, dp=NULL)
}


\arguments{
\item{x}{vector of quantiles. Missing values (\code{NA}s) are allowed.}
\item{p}{vector of probabililities.}
\item{xi}{vector of location parameters.}
\item{omega}{vector of scale parameters; must be positive.}
\item{alpha}{vector of slant parameters. With \code{pst} and \code{qst}, 
it must be of length 1.}
\item{nu}{a single positive value representing the degrees of freedom;
  it can be non-integer. Default value is \code{nu=Inf} which corresponds 
  to the skew-normal distribution.
}

\item{dp}{a vector of length 4, whose elements represent location, scale
  (positive), slant and degrees of freedom, respectively.  If \code{dp} is
  specified, the individual parameters cannot be set.}

\item{n}{a positive integer representing the sample size.}

\item{log, log.p}{logical; if \code{TRUE}, densities  are given as log-densities
   and probabilities \code{p} are given as \code{log(p)}}

\item{tol}{
  a scalar value which regulates the accuracy of the result of
  \code{qsn}, measured on the probability scale.
  }

\item{method}{an integer value between \code{0} and \code{5} which selects 
  the computing method; see \sQuote{Details} below for the meaning of these
  values. If \code{method=0} (default value), an automatic choice is made
  among the four actual computing methods, depending on the other
  arguments.}

\item{lower.tail}{
  logical; if \code{TRUE} (default), probabilities are \eqn{P\{X\le x\}},
  otherwise  \eqn{P\{X\ge x\}}
  }

\item{...}{additional parameters passed to \code{integrate}; 
           see \sQuote{Details} for further information.}

}

\value{Density (\code{dst}), probability (\code{pst}), quantiles (\code{qst}) 
and random sample (\code{rst}) from the skew-\eqn{t} distribution with given 
\code{xi}, \code{omega}, \code{alpha} and \code{nu} parameters.}

\section{Details}{
Typical usages are
\preformatted{%
dst(x, xi=0, omega=1, alpha=0, nu=Inf, log=FALSE)
dst(x, dp=, log=FALSE)
pst(x, xi=0, omega=1, alpha=0, nu=Inf, method=0, ...)
pst(x, dp=, log=FALSE)
qst(p, xi=0, omega=1, alpha=0, nu=Inf, tol=1e-8,  method=0, ...)
qst(x, dp=, log=FALSE)
rst(n=1, xi=0, omega=1, alpha=0, nu=Inf)
rst(x, dp=, log=FALSE)
}
}
\section{Background}{
The family of skew-\eqn{t} distributions is an extension of the Student's
\eqn{t} family, via the introduction of a \code{alpha} parameter which 
regulates skewness; when \code{alpha=0}, the skew-\eqn{t} distribution 
reduces to the usual Student's \eqn{t} distribution. 
When \code{nu=Inf}, it reduces to the skew-normal distribution. 
When \code{nu=1}, it reduces to a form of skew-Cauchy distribution.
See Chapter 4 of Azzalini & Capitanio (2014) for additional information. 
A multivariate version of the distribution exists; see \code{dmst}.
}

\section{Details}{
For evaluation of \code{pst}, and so indirectly of
\code{qst}, four different methods are employed.
In all the cases, the actual computations are performed for the normalized 
values \code{z=(x-xi)/omega}).
Method 1 consists in using \code{pmst} with dimension \code{d=1}.
Method 2 applies \code{integrate} to the density function \code{dst}.
Method 3 again uses \code{integrate} too but with a different integrand,
as given in Section 4.2 of Azzalini & Capitanio (2003, full version of
the paper).
Method 4 consists in the recursive procedure of Jamalizadeh, Khosravi and
Balakrishnan (2009), which is recalled in Complement 4.3 on 
Azzalini & Capitanio (2014); the recursion over \code{nu} starts from 
the explicit expression for \code{nu=1} given by \code{psc}.
Method 5 is targeted to tail probabilities only, and it returns \code{NA}s
for non-extreme \code{x} values (those with \code{abs(z)<=20}); 
it is based on expressions given in Complement 4.4 of Azzalini and Capitanio (2014).
Method 1 and 4 are only suitable for integer values of \code{nu}.
Method 4 becomes progressively less efficient as \code{nu} increases,
because the value of \code{nu} determines the number of nested calls, 
but the decay of efficiency is slower for larger values of \code{length(x)}.
If the default argument value \code{method=0} is retained, an automatic choice
among the above four methods is made, which depends on the values of 
\code{nu, alpha, z}.  

The \code{pst} argument \code{...} can be used to set values of \code{rel.tol} 
and \code{rel.tol} passed to \code{integrate}, which is used if \code{method=2}
or \code{method=3}, or possibly when \code{method=0} selects one of these options.
For other values of \code{method}, use of \code{...} has no effect. 
This mechanism also regulates \code{qst}, since this function hinges on \code{pst}.

If \code{qst} is called with \code{nu>1e4}, the computation is transferred to
\code{qsn}. 
}

\references{
Azzalini, A. and Capitanio, A. (2003).
 Distributions generated by perturbation of symmetry 
 with emphasis on a multivariate skew-\emph{t} distribution.
 \emph{J.Roy. Statist. Soc. B} \bold{65}, 367--389.
 Full version of the paper at \url{https://arXiv.org/abs/0911.2342}.

Azzalini, A. with the collaboration of Capitanio, A. (2014). 
 \emph{The Skew-normal and Related Families}. 
 Cambridge University Press, IMS Monographs series. 

Jamalizadeh, A., Khosravi, M., and Balakrishnan, N. (2009).
  Recurrence relations for distributions of a skew-\emph{t} and a linear
  combination of order statistics from a bivariate-\emph{t}.
  \emph{Comp. Statist. Data An.} \bold{53}, 847--852.
}
  

\seealso{\code{\link{dmst}}, \code{\link{dsn}}, \code{\link{dsc}}, 
         \code{\link[stats]{integrate}}}

\examples{
pdf <- dst(seq(-4, 4, by=0.1), alpha=3, nu=5)
rnd <- rst(100, 5, 2, -5, 8)
q <- qst(c(0.25, 0.50, 0.75), alpha=3, nu=5)
pst(q, alpha=3, nu=5)  # must give back c(0.25, 0.50, 0.75)
#
p1 <- pst(x=seq(-3,3, by=1), dp=c(0,1,pi, 3.5))
p2 <- pst(x=seq(-3,3, by=1), dp=c(0,1,pi, 3.5), method=2, rel.tol=1e-9)
}
\keyword{distribution}

