\name{genhypergeo}
\alias{genhypergeo}
\alias{genhypergeo_series}
\alias{genhypergeo_contfrac}
\title{The generalized hypergeometric function}
\description{
  The generalized hypergeometric function, using either the series
  expansion or the continued fraction expansion.
}
\usage{
genhypergeo(U, L, z, tol=0, maxiter=2000, check_mod=TRUE,
    polynomial=FALSE, debug=FALSE, series=TRUE)
genhypergeo_series(U, L, z, tol=0, maxiter=2000, check_mod=TRUE,
    polynomial=FALSE, debug=FALSE) 
genhypergeo_contfrac(U, L, z, tol = 0, maxiter = 2000)
}
\arguments{
  \item{U,L}{Upper and lower arguments respectively (real or complex)}
  \item{z}{Primary complex argument (see notes)}
  \item{tol}{tolerance with default zero meaning to iterate until
    additional terms to not change the partial sum}
  \item{maxiter}{Maximum number of iterations to perform}
  \item{check_mod}{Boolean, with default \code{TRUE} meaning to check
    that the modulus of \code{z} is less than 1}
  \item{polynomial}{Boolean, with default \code{FALSE} meaning to
    evaluate the series until converged, or return a warning; and
    \code{TRUE} meaning to return the sum of \code{maxiter} terms,
    whether or not converged.  This is useful when either \code{A}
    or \code{B} is a nonpositive integer in which case the hypergeometric
    function is a polynomial}
  \item{debug}{Boolean, with \code{TRUE} meaning to return debugging
    information and default \code{FALSE} meaning to return just the
    evaluate} 
  \item{series}{In function \code{genhypergeo()}, Boolean argument with
    default \code{TRUE} meaning to return the result of
    \code{genhypergeo_series()} and \code{FALSE} the result of
    \code{genhypergeo_contfrac()}} 
}
\details{
  Function \code{genhypergeo()} is a wrapper for functions
  \code{genhypergeo_series()} and \code{genhypergeo_contfrac()}.
  
  Function \code{genhypergeo_series()} is the workhorse for the whole
  package; every call to \code{hypergeo()} uses this function except for
  the (apparently rare---but see the examples section) cases where
  continued fractions are used.

  The generalized hypergeometric function [here \code{genhypergeo()}]
  appears from time to time in the literature (eg Mathematica) as
  
  \deqn{F(U,L;z) = \sum_{n=0}^\infty\frac{(u_1)_n(u_2)_n\ldots
      (u_i)_n}{(l_1)_n(l_2)_n\ldots
      (l_j)_n}\cdot\frac{z^n}{n!}}{[omitted; see PDF]} where
      \eqn{U=\left(u_1,\ldots,u_i\right)}{U=(u_1,...,u_i)} and
      \eqn{L=\left(l_1,\ldots,l_i\right)}{L=(l_1,...,l_i)} are the
      \dQuote{upper} and \dQuote{lower} vectors respectively.  The
      radius of convergence of this formula is 1.
  
  For the Confluent Hypergeometric function, use \code{genhypergeo()} with
  length-1 vectors for arguments \code{U} and \code{V}.
  
  For the \eqn{{}_0\!F_1}{0F1} function (ie no \dQuote{upper} arguments), use
  \code{genhypergeo(NULL,L,x)}.
    
  See documentation for \code{genhypergeo_contfrac()} for details of
  the continued fraction representation.
}
\references{
  M. Abramowitz and I. A. Stegun 1965. \emph{Handbook of
    mathematical functions}.  New York: Dover
}
\author{Robin K. S. Hankin}
\note{
  The radius of convergence for the series is 1 but under some
  circumstances, analytic continuation defines a function over the whole
  complex plane (possibly cut along \eqn{(1,\infty)}{(1,inf)}).  Further
  work would be required to implement this.
}
\seealso{\code{\link{hypergeo}},\code{\link{genhypergeo_contfrac}}}
\examples{

genhypergeo(U=c(1.1,0.2,0.3), L=c(10.1,pi*4), check_mod=FALSE, z=1.12+0.2i)
genhypergeo(U=c(1.1,0.2,0.3), L=c(10.1,pi*4),z=4.12+0.2i,series=FALSE)

}
\keyword{math}
