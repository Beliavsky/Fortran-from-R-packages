\name{hypergeo_cover1}
\alias{hypergeo_cover1}
\alias{hypergeo_cover2}
\alias{hypergeo_cover3}
\title{Hypergeometric functions for special values of the parameters}
\description{
Hypergeometric functions for special values of the parameters
}
\usage{
hypergeo_cover1(A, B, m, z, tol = 0, maxiter = 2000, method = "a", give = FALSE)
hypergeo_cover2(A, C, m, z, tol = 0, maxiter = 2000, method = "a", give = FALSE)
hypergeo_cover3(A, n, m, z, tol = 0, maxiter = 2000, method = "a", give = FALSE)
}
\arguments{
  \item{A,B,C}{parameters for the hypergeometric function}
  \item{m,n}{Integers (positive or negative)}
  \item{z}{Primary complex argument}
  \item{tol,maxiter}{Numerical arguments passed to \code{genhypergeo()}}
  \item{method}{Method, passed to \code{f15.3.10()} (qv)}
  \item{give}{Boolean with \code{TRUE} meaning to return the choice of
    helper function used (eg \code{f15.3.7()}), and default \code{FALSE}
    meaning to return the hypergeometric function's value}
}
\details{

  These functions deal with the exceptional cases listed on page
  559-560.

  \itemize{
    \item Function \code{hypergeo_cover1()} deals with the case
    \eqn{C=A+B\pm m,m=0,1,2,\ldots}{C=A+B +/- m, m=0,1,2,...}
    \item Function \code{hypergeo_cover2()} deals with the case
    \eqn{B=A\pm m,m=0,1,2,\ldots}{B=A +/- m, m=0,1,2,...}
    \item Function \code{hypergeo_cover3()} deals with the case
    \eqn{C-A=0,-1,-2,\ldots}{C-A=0,-1,-2,...} [elementary] and
    \eqn{C-A=1,2,\ldots}{C-A=1,2,...} [not covered by AMS-55]
    }
}
\references{
  M. Abramowitz and I. A. Stegun 1965. \emph{Handbook of
    mathematical functions}.  New York: Dover
}
\note{
  Function \code{hypergeo_cover3()} is required because the
    \dQuote{limiting process} mentioned on p560, just after 15.3.14, is
    not explicit.  Which is why it dispatches to
    \code{w07.23.06.0026.01()} and \code{w07.23.06.0031.01()},
    documented at \code{wolfram}.
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{hypergeo}},\code{\link{f15.3.10}},\code{\link{wolfram}}}
\examples{


# Test hypergeo_cover1():
jjR <- hypergeo(pi,pi/2,3*pi/2-4, z=0.1+0.2i)
jjM <- 0.53745229690249593045 + 1.8917456473240515664i


# Test hypergeo_cover2():
jjM <- -0.15888831928748121465e-5 + 0.40339599711492215912e-4i
jjR <- hypergeo(pi,pi+2, 1.1 , 1+10i)  # This is 15.3.13
stopifnot(Mod(jjR-jjM)<1e-10)


# Test hypergeo_cover3()
jjM <- -0.24397135980533720308e-1 + 0.28819643319432922231i
jjR <- hypergeo(pi, 1.4, pi+4, 1+6i)
stopifnot(Mod(jjR-jjM)<1e-10)

}
\keyword{math}
