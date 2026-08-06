%  file sn/man/zeta.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 1998,2013 Adelchi Azzalini
%---------------------
\name{zeta}
\alias{zeta}


\title{Function \eqn{\log(2\,\Phi(x))}{log(2 Phi(x))} and its derivatives}

\description{The function \code{log(2*pnorm(x))} and its derivatives, 
  including the inverse Mills ratio.}
  
\usage{zeta(k, x)}

\arguments{
\item{k}{an integer number between 0 and 5.}

\item{x}{a numeric vector. Missing values (\code{NA}s)  and \code{Inf}s are 
  allowed.}
}

\value{a vector representing the \code{k}-th order derivative evaluated at 
  \code{x}.}

\details{
  For \code{k} between 0 and 5, the derivative of  order \code{k}  of
  \eqn{\log(2\,\Phi(x))}{log(2\Phi(x))} is  evaluated, where \eqn{\Phi(x)} 
  denotes the \eqn{N(0,1)} cumulative distribution function.
  The derivative of order \code{k=0} refers to the function itself.
  If \code{k} is not an integer within \code{0,..., 5}, \code{NULL} is returned.
}

\section{Background}{
  The computation for \code{k>1} is reduced to the case \code{k=1}, making use
  of expressions given by Azzalini and Capitanio (1999); see especially 
  Section 4 of the full-length version of the paper.   The main facts are 
  summarized in Section 2.1.4 of Azzalini and Capitanio (2014).
  
  For numerical stability, the evaluation of \code{zeta(1,x)} when 
  \code{x < -50} makes use of the asymptotic expansion (26.2.13) of 
  Abramowitz and Stegun (1964).

  \code{zeta(1,-x)} equals \code{dnorm(x)/pnorm(-x)} (in principle, apart from
  the above-mentioned asymptotic expansion), called the  
  \emph{inverse Mills ratio}.
}

\references{
Abramowitz, M. and Stegun, I. A., editors (1964).
  \emph{Handbook of Mathematical Functions}. 
  Dover Publications.

Azzalini, A. and Capitanio, A. (1999).
  Statistical applications of the multivariate skew normal distribution.
  \emph{J.Roy.Statist.Soc. B} \bold{61}, 579--602. Full-length version 
  available at \url{https://arXiv.org/abs/0911.2093}
     
Azzalini, A. with the collaboration of Capitanio, A. (2014). 
 \emph{The Skew-Normal and Related Families}. 
 Cambridge University Press, IMS Monographs series. 
}
\examples{
y <- zeta(2,seq(-20,20,by=0.5))
#
for(k in 0:5) curve(zeta(k,x), from=-1.5, to=5, col = k+2, add = k > 0)
legend(3.5, -0.5, legend=as.character(0:5), col=2:7, lty=1)
}

\keyword{math}
\concept{Mills ratio}
