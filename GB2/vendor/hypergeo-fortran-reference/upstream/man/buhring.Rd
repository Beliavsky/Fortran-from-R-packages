\name{buhring}
\alias{buhring}
\alias{buhring_eqn11}
\alias{buhring_eqn12}
\alias{buhring_eqn5_factors}
\alias{buhring_eqn5_series}
\alias{hypergeo_buhring}
\alias{buhring_eqn12}
\title{Evaluation of the hypergeometric function using Buhring's method}
\description{
Expansion of the hypergeometric function
using the residue theorem; useful for when the primary argument is close
to the critical points \eqn{1/2\pm i\sqrt{3}/2}{0.5+/-i.sqrt(3)/2}
}
\usage{
hypergeo_buhring(A,B,C,z,z0=1/2,tol=0,maxiter=2000,use11=TRUE)
buhring_eqn11(n,S,A,B,C,z0=1/2)
buhring_eqn12(n,S,A,B,C,z0=1/2)
buhring_eqn5_factors(A,B,C,z,z0=1/2)
buhring_eqn5_series(S,A,B,C,z,z0=1/2,use11=FALSE,tol=0,maxiter=2000)
}
\arguments{
  \item{A,B,C}{Parameters (real)}
  \item{S}{Parameter taken to be either \code{A} or \code{B}}
  \item{n}{Term to calculate in \code{buhring_eqn11()} or
    \code{buhring_eqn12()}}
  \item{z}{Primary complex argument}
  \item{z0}{Centre of circle of non-convergence; series expressed in
    powers of \eqn{1/\left(z-z_0\right)^n}{1/(z-z0)^n}}
  \item{tol,maxiter}{tolerance and maximum number of iterations (as
    in \code{hypergeo()})}
  \item{use11}{Boolean with default \code{TRUE} meaning to use
    \code{buhring_eqn11()} and \code{FALSE} meaning to use \code{buhring_eqn12()}}
}
\details{

The functions are direct transcriptions of Buhring 1987.  The basic
idea is to expand the hypergeometric function in powers of
\eqn{(z-z_0)^{-1}}{1/(z-z0)}.

Functions \code{buhring_eqn11()} and \code{buhring_eqn12()} return the
coefficients \eqn{d_n}{d(n)} given by equations 11 and 12 of Buhring
1987.  

The series do not converge satisfactorily near the critical points due
to some sort of numerical instability.  But they seem to work OK if
\eqn{|z-1/2|} is large.

}
\references{ 
  \itemize{
    \item W. Buhring 1987.  \dQuote{An analytic continuation of the
      hypergeometric series}, \emph{Siam J. Math. Anal.} 18(3)
  }
}
\author{Robin K. S. Hankin}
\note{
There is some issue which prevents the series from converging correctly,
also sometimes the sequence converges to a demponstrably incorrect
value. 
}
\seealso{\code{\link{residue}}}
\examples{


# should be identical:
buhring_eqn11(n=0:10,S=1/2,A=1/2,B=1/3,C=pi)
buhring_eqn12(n=0:10,S=1/2,A=1/2,B=1/3,C=pi)
# but differ in one element


a <- hypergeo(1/2,1/3,4,1+8i,maxiter=90)
b <- hypergeo_buhring(1/2,1/3,4,1+8i,maxiter=90)
# should be identical but are not 


# following command fails due to numerical instability:
\dontrun{
hypergeo_buhring(1/2,1/3,pi,z=1/2 + 1i*sqrt(3)/2)
}

}
\keyword{math}

