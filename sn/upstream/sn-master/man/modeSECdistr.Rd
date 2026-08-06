%  file sn/man/mmodeSECdistr.Rd  
%  This file is a component of the package 'sn' for R
%  copyright (C) 2013 Adelchi Azzalini
%---------------------
\name{modeSECdistr}
\alias{modeSECdistr}
\title{The mode of a skew-elliptically contoured (\acronym{SEC}) distribution}

\description{
  Compute compute the mode of a univariate or multivariate \acronym{SEC}
 distribution.}

\usage{
  modeSECdistr(dp, family, object=NULL)
}

\arguments{

\item{dp}{
   a numeric vector (in the univariate case, for class \code{SECdistrUv}) 
   or a list (in the multivariate case, , for class \code{SECdistrUv}) 
   of parameters which identify the specific distribution within the named
   \code{family}.
   }
    
\item{family}{a character string which identifies the parametric
   family among those admissible for classes \code{SECdistrUv} or 
   \code{SECdistrMv}.}

\item{object}{an object of class \code{SECdistrUv} or \code{SECdistrMv} as
   created by \code{\link{makeSECdistr}} or \code{\link{extractSECdistr}};
   if this argument is used, arguments \code{dp} and \code{family} must 
   not be set, and \emph{vice versa}.
   }
}  

      
\value{a numeric vector}

\section{Background}{The mode is obtained through numerical maximization.
In the multivariate case, the problem is reduced to a one-dimensional search 
using Propositions 5.14 and 6.2 of the reference below.}

\references{
 Azzalini, A. with the collaboration of Capitanio, A. (2014). 
 \emph{The Skew-Normal and Related Families}. 
 Cambridge University Press, IMS Monographs series. 
}  
 
\seealso{\code{\link{makeSECdistr}} and \code{\link{extractSECdistr}}
  for additional information and for constructing a suitable \code{object},
   
 \code{\link{SECdistrUv-class}} and \code{\link{SECdistrMv-class}}
 for methods \code{mean} and \code{vcov} which compute the mean (vector) 
 and the variance (matrix) of the  \code{object} distribution
 }
 

  
\examples{
dp3 <- list(xi=1:3, Omega=toeplitz(1/(1:3)), alpha=c(3,-1,2), nu=5)
st3 <- makeSECdistr(dp3, family="ST", name="ST3", compNames=c("U", "V", "W"))
mode1 <- modeSECdistr(dp3, "ST")
mode2 <- modeSECdistr(object=st3) # the same of mode1
}

\keyword{multivariate}
\keyword{distribution}
