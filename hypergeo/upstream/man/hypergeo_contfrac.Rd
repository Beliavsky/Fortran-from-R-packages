\name{hypergeo_contfrac}
\alias{hypergeo_contfrac}
\alias{genhypergeo_contfrac_single}
\title{Continued fraction expansion of the hypergeometric function}
\description{
Continued fraction expansion of the hypergeometric and generalized
hypergeometric functions using continued fraction expansion.
}
\usage{
hypergeo_contfrac(A, B, C, z, tol = 0, maxiter = 2000)
genhypergeo_contfrac_single(U, L, z, tol = 0, maxiter = 2000)
}
\arguments{
  \item{A,B,C}{Parameters (real or complex)}
  \item{U,L}{In function \code{genhypergeo_contfrac()}, upper and lower
    arguments as in \code{genhypergeo()}} 
  \item{z}{Complex argument}
  \item{tol}{tolerance (passed to \code{GCF()})}
  \item{maxiter}{maximum number of iterations}
}
\details{
  These functions are included in the package in the interests of
  completeness, but it is not clear when it is advantageous to use
  continued fraction form rather than the series form.
}
\references{
  \itemize{
    \item
    M. Abramowitz and I. A. Stegun 1965. \emph{Handbook of
      mathematical functions}.  New York: Dover
    \item \code{http://functions.wolfram.com/Hypergeometric2F1.pdf}
  }
}
\author{Robin K. S. Hankin}
\note{

  The continued fraction expression is the RHS identity 07.23.10.0001.01 of
  \code{Hypergeometric2F1.pdf}.
  
  The function sometimes fails to converge to the correct value but no
  warning is given.

  Function \code{genhypergeo_contfrac()} is documented under
  \code{genhypergeo.Rd}. 
}
\seealso{\code{\link{genhypergeo}}}
\examples{
hypergeo_contfrac(0.3 , 0.6 , 3.3 , 0.1+0.3i)
# Compare Maple: 1.0042808294775511972+0.17044041575976110947e-1i

genhypergeo_contfrac_single(U=0.2 , L=c(9.9,2.7,8.7) , z=1+10i)
# (powerseries does not converge)
# Compare Maple: 1.0007289707983569879 + 0.86250714217251837317e-2i

}
\keyword{math}
