%  file sn/man/Qpenalty.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{Qpenalty}
\alias{Qpenalty}
\alias{MPpenalty}
\concept{penalized likelihood}
\concept{prior distribution}

\title{Penalty function for log-likelihood of \code{selm} models}

\description{Penalty function for the log-likelihood of \code{selm} models 
  when \code{method="MPLE"}. \code{Qpenalty} is the default function;
  \code{MPpenalty} is an example of a user-defined function effectively
  corresponding to a prior distributio on \code{alpha}. }

\usage{
Qpenalty(alpha_etc, nu = NULL, der = 0)

MPpenalty(alpha, der = 0)    
}

\arguments{

 \item{alpha_etc, alpha}{in the univariate case, a single value \code{alpha}; 
   in the multivariate case, a two-component list whose first component is 
   the vector \code{alpha}, the second one is matrix \code{cov2cor(Omega)}. }

 \item{nu}{degrees of freedom, only required if \code{selm} is called
   with \code{family="ST"}. }

  \item{der}{a numeric value in the set \kbd{0,1,2} which indicates the
   required numer of derivatives of the function. In the multivariate case
   the function will only be called with \code{der} equal to 0 or 1.}
}

\details{The penalty is a function of \code{alpha}, but its expression may
depend on other ingredients, specifically \code{nu} and \code{cov2cor(Omega)}.
See \sQuote{Details} of \code{\link{selm}} for additional information.

The penalty mechanism allows to introduce a prior distribution \eqn{\pi} 
for \eqn{\alpha} by setting \eqn{Q=-\log\pi}{Q=-log(\pi)}, 
leading to a maximum \emph{a posteriori} estimate in the stated sense. 

As a simple illustration of this mechanism, function \code{MPpenalty}
implements the `matching prior' distribution for the univariate \acronym{SN}
distribution studied by Cabras \emph{et al.} (2012); a brief summary of the
proposal is provided in Section 3.2 of Azzalini and Capitanio (2014). Note
that, besides \code{alpha=+/-Inf}, this choice also penalizes \code{alpha=0}
with \code{Q=Inf}, effectively removing \code{alpha=0} from the parameter
space.  

Starting from the code of function \code{MPpenalty}, a user should be able 
to introduce an alternative prior distribution if so desired.
}

\value{A positive number \code{Q} representing the penalty, possibly
  with attributes \code{attr(Q, "der1")} and \code{attr(Q, "der2")}, 
  depending onthe input value \code{der}.} 


\references{
Azzalini, A. with the collaboration of Capitanio, A. (2014). 
  \emph{The Skew-Normal and Related Families}. 
  Cambridge University Press, IMS Monographs series.
 
Cabras, S., Racugno, W., Castellanos, M. E., and Ventura, L. (2012).
  A matching prior for the shape parameter of the skew-normal distribution.
  \emph{Scand. J. Statist.} \bold{39}, 236--247.
} 

\author{Adelchi Azzalini}

% \note{%%  ~~further notes~~}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{selm}} function}


\examples{
data(frontier)
m2 <- selm(frontier ~ 1)  # no penalty
m2a <- selm(frontier ~ 1, method="MPLE") # penalty="Qpenalty" is implied here
m2b <- selm(frontier ~ 1, method="MPLE", penalty="MPpenalty")    
}
 

%\keyword{ ~kwd1 }
