%  file sn/man/dsc.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{dsc}
\alias{dsc}
\alias{psc}
\alias{qsc}
\alias{rsc}
\title{Skew-Cauchy Distribution}

\description{Density function, distribution function, quantiles and random
 number generation for the skew-Cauchy (SC) distribution.}

\usage{
dsc(x, xi = 0, omega = 1, alpha = 0, dp = NULL, log = FALSE)
psc(x, xi = 0, omega = 1, alpha = 0, dp = NULL)
qsc(p, xi = 0, omega = 1, alpha = 0, dp = NULL) 
rsc(n = 1, xi = 0, omega = 1, alpha = 0, dp = NULL)
}

\arguments{
 \item{x}{vector of quantiles. Missing values (\code{NA}s) and \code{Inf}'s
   are allowed.}

 \item{p}{vector of probabilities. Missing values (\code{NA}s) are allowed.}

 \item{xi}{ vector of location parameters.}

 \item{omega}{vector of (positive) scale parameters.}

 \item{alpha}{vector of slant parameters.}

 \item{dp}{a vector of length 3 whose elements represent the parameters 
   described above. If \code{dp} is specified, the individual parameters 
   cannot be set.}

 \item{n}{sample size.}

 \item{log}{logical flag used in \code{dsc} (default \code{FALSE}).
   When \code{TRUE}, the logarithm of the density values is returned.}

}

\value{density (\code{dsc}), probability (\code{psc}), quantile (\code{qsc})
  or random sample (\code{rsc}) from the skew-Cauchy distribution with given
  \code{xi}, \code{omega} and \code{alpha} parameters or from the extended
  skew-normal if \code{tau!=0} }

\section{Details}{
Typical usages are
\preformatted{%
dsc(x, xi=0, omega=1, alpha=0, log=FALSE)
dsc(x, dp=, log=FALSE)
psc(x, xi=0, omega=1, alpha=0)
psc(x, dp= )
qsc(p, xi=0, omega=1, alpha=0)
qsc(x, dp=)
rsc(n=1, xi=0, omega=1, alpha=0)
rsc(x, dp=)
}
}

\section{Background}{
The skew-Cauchy distribution can be thought as a skew-\eqn{t} with tail-weight
parameter \code{nu=1}. In this case, closed-form expressions of the 
distribution function and the quantile function have been obtained by
Behboodian \emph{et al.} (2006).
The key facts are summarized in Complement 4.2 of Azzalini and Capitanio (2014).
A multivariate version of the distribution exists.
}

\references{
 
Azzalini, A. with the collaboration of Capitanio, A. (2014). 
\emph{The Skew-normal and Related Families}. 
Cambridge University Press, IMS Monographs series.

Behboodian, J., Jamalizadeh, A., and Balakrishnan, N. (2006).
A new class of skew-Cauchy distributions.
\emph{Statist. Probab. Lett.} \bold{76}, 1488--1493.

}

\seealso{\code{\link{dst}},  \code{\link{dmsc}}}

\examples{
pdf <- dsc(seq(-5,5,by=0.1), alpha=3)
cdf <- psc(seq(-5,5,by=0.1), alpha=3)
q <- qsc(seq(0.1,0.9,by=0.1), alpha=-2)
p <- psc(q, alpha=-2)  
rn <- rsc(100, 5, 2, 5)
}
\keyword{distribution}
