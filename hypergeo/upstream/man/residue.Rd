\name{residue}
\alias{residue}
\alias{hypergeo_residue}
\alias{hypergeo_residue_general}
\alias{hypergeo_residue_close_to_crit_single}
\alias{hypergeo_residue_close_to_crit_multiple}

\title{Evaluation of the hypergeometric function using the residue theorem}
\description{
Expansion of the hypergeometric function
using the residue theorem; useful for when the primary argument is close
to the critical points \eqn{1/2\pm i\sqrt{3}/2}{0.5+/-i.sqrt(3)/2}
}
\usage{
hypergeo_residue_general(A, B, C, z, r, O=z, tol=0, maxiter=2000)
hypergeo_residue_close_to_crit_single(A, B, C, z, strategy='A', tol=0, maxiter=2000)
hypergeo_residue_close_to_crit_multiple(A, B, C, z, strategy='A', tol=0, maxiter=2000)
}
\arguments{
  \item{A,B,C}{Parameters (real or complex)}
  \item{z}{Complex argument}
  \item{tol,maxiter}{tolerance and maximum number of iterations (passed to \code{hypergeo()})}
  \item{r,O}{Radius and center of circle to integrate over}
  \item{strategy}{Indicates which strategy to use.  Strategy \sQuote{A}
    means to use the critical point as the centre of the circle and
    strategy \sQuote{B} means to use \eqn{z}}
}
\details{

These functions are not really intended for the user; \code{hypergeo()}
uses \code{hypergeo_residue_close_to_crit_multiple()} when
\eqn{\left|z-c\right|}{abs(z-c)} is less than \eqn{0.1} (hardwired) for
\eqn{c} being either of the two critical points.  Infinite regress is
avoided because the contour is always more than this distance from the
critical points.
  
These functions use the residue theorem
\eqn{f\left(z_0\right)=\oint_C\frac{f(z)\,dz}{z-z_0}}{f(z_0)=int
f(z)dz/(z-z_0)} to evaluate the hypergeometric function near the two
critical points \eqn{1/2\pm i\sqrt{3}/2}{0.5+/-i.sqrt(3)/2}.  These
points are problematic because all of the transformations listed under
\code{thingfun()} take the points either to themselves or each other.

At these points the ratio of successive terms in the hypergeometric
series tends to one and thus numerical summation is difficult.

The hypergeometric function, however, is not at all badly behaved near
these critical points (see examples); but OTOH there do not seem to be
any identities for the hypergeometric function at these points.

I have not investigated in detail whether strategy \sQuote{A} or
\sQuote{B} is better.  I would expect that \sQuote{A} is faster but
\sQuote{B} more accurate, on the grounds that \sQuote{A} uses a contour
whose closest approach to the critical point is further than that of
\sQuote{B}; but \sQuote{B} uses a contour which does not vary in
distance from \eqn{z}.

But both seem to be fairly accurate and fairly fast, and I have not
systematically investigated the pros and cons.
}
\references{ 
  \itemize{
    \item W. Buhring 1987.  \dQuote{An analytic continuation of the
      hypergeometric series}, \emph{Siam J. Math. Anal.} 18(3)
  }
}
\author{Robin K. S. Hankin}
\note{
The residue theorem appears to be absurdly accurate for numerical evaluation
}
\seealso{\code{\link{buhring}}}
\examples{
c1 <- 1/2-sqrt(3)/2i
c2 <- 1/2+sqrt(3)/2i

a1_R <- hypergeo(1/2,1/3,pi,c1)
a1_M <- 1.0154051314906669 + 0.0544835896509068i

x <- y <- seq(from=-0.1,to=0.1,len=100)
elliptic::view(x,y,hypergeo(1/2,1,1/3,outer(x,1i*y,"+")))

}
\keyword{math}

