\name{hypergeo}
\alias{hypergeo}
\title{The hypergeometric function}
\description{
The Hypergeometric and generalized hypergeometric functions as defined
by Abramowitz and Stegun.    Function \code{hypergeo()} is the user
interface to the majority of the package functionality; it dispatches to
one of a number of subsidiary functions.
}
\usage{
hypergeo(A, B, C, z, tol = 0, maxiter=2000)
}
\arguments{
  \item{A,B,C}{Parameters for \code{hypergeo()}}
  \item{z}{Primary argument, complex}
  \item{tol}{absolute tolerance; default value of zero means to continue
  iterating until the result does not change to machine precision;
  strictly positive values give less accuracy but faster evaluation}
  \item{maxiter}{Integer specifying maximum number of iterations}
}
\details{
  The hypergeometric function as defined by Abramowitz and Stegun,
  equation 15.1.1, page 556 is
  \deqn{ {}_2F_1(a,b;c;z) =
    \sum_{n=0}^\infty\frac{(a)_n(b)_n}{(c)_n}\cdot\frac{z^n}{n!}}{[omitted;
    see PDF]}
  
  where
  \eqn{(a)_n=a(a+1)\ldots(a+n-1)=\Gamma(a+n)/\Gamma(a)}{(a)_n=Gamma(a+n)/Gamma(a)}
  is the Pochammer symbol (6.1.22, page 256).
  
  Function \code{hypergeo()} is the front-end for a rather unwieldy set
  of back-end functions which are called when the parameters \code{A},
  \code{B}, \code{C} take certain values.

  The general case (that is, when the parameters do not fall into a
  \dQuote{special} category), is handled by \code{hypergeo_general()}.
  This applies whichever of the transformations given on page 559 gives
  the smallest modulus for the argument \code{z}.

  Sometimes \code{hypergeo_general()} and all the transformations on
  page 559 fail to converge, in which case \code{hypergeo()} uses the
  continued fraction expansion \code{hypergeo_contfrac()}.
  
  If this fails, the function uses integration via \code{f15.3.1()}.
  
}
\references{Abramowitz and Stegun 1955.
  \emph{Handbook of mathematical functions with formulas, 
      graphs and mathematical tables} (AMS-55).
      National Bureau of Standards}
\author{Robin K. S. Hankin}
\note{
Abramowitz and Stegun state:

\dQuote{The radius of convergence of the Gauss hypergeometric series
\eqn{\ldots}{...} is \eqn{\left|z\right|=1}{|z|=1}} (AMS-55, section
15.1, page 556).

This reference book gives the correct radius of convergence; use the
ratio test to verify it.  Thus if \eqn{|z|>1}, the hypergeometric series
will diverge and function \code{genhypergeo()} will fail to converge.

However, the hypergeometric function is defined over the whole of the
complex plane, so analytic continuation may be used if appropriate cut
lines are used.  A cut line must join \eqn{z=1} to (complex) infinity;
it is conventional for it to follow the real axis in a positive
direction from \eqn{z=1} but other choices are possible.

Note that in using the package one sometimes draws a
\dQuote{full precision not achieved} warning from \code{gamma()}; and
complex arguments are not allowed.  I would suggest either ignoring the
warning (the error of \code{gamma()} is unlikely to be large) or to use
one of the bespoke functions such as \code{f15.3.4()} and tolerate the
slower convergence, although this is not always possible.
}
\seealso{\code{\link{hypergeo_powerseries}},
\code{\link{hypergeo_contfrac}}, \code{\link{genhypergeo}}}

\examples{
#  equation 15.1.3, page 556:
f1 <- function(x){-log(1-x)/x}
f2 <- function(x){hypergeo(1,1,2,x)}
f3 <- function(x){hypergeo(1,1,2,x,tol=1e-10)}
x <- seq(from = -0.6,to=0.6,len=14)
f1(x)-f2(x)
f1(x)-f3(x)  # Note tighter tolerance

# equation 15.1.7, p556:
g1 <- function(x){log(x + sqrt(1+x^2))/x}
g2 <- function(x){hypergeo(1/2,1/2,3/2,-x^2)}
g1(x)-g2(x)  # should be small 
abs(g1(x+0.1i) - g2(x+0.1i))  # should have small modulus.

# Just a random call, verified by Maple [ Hypergeom([],[1.22],0.9087) ]:
genhypergeo(NULL,1.22,0.9087)


# Little test of vectorization (warning: inefficient):
hypergeo(A=1.2+matrix(1:10,2,5)/10, B=1.4, C=1.665, z=1+2i)


# following calls test for former bugs:
hypergeo(1,2.1,4.1,1+0.1i)
hypergeo(1.1,5,2.1,1+0.1i)
hypergeo(1.9, 2.9, 1.9+2.9+4,1+0.99i) # c=a+b+4; hypergeo_cover1()


}
\keyword{math}

