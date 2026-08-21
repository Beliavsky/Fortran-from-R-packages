\name{is.nonpos}
\alias{is.nonpos}
\alias{is.near_integer}
\alias{is.zero}
\alias{isgood}
\alias{thingfun}
\alias{crit}
\alias{lpham}
\title{Various utilities}
\description{
Various utilities needing nonce functions
}
\usage{
is.near_integer(i, tol=getOption("tolerance"))
is.nonpos(i)
is.zero(i)
isgood(x, tol)
thingfun(z, complex=FALSE)
crit(...)
lpham(x,n)
}
\arguments{
  \item{i}{Numerical vector of suspected integers}
  \item{tol}{Tolerance}
  \item{x}{Argument to \code{isgood()} and \code{lpham()}}
  \item{z}{Complex vector}
  \item{complex}{In function \code{thingfun()}, Boolean with default
    \code{FALSE} meaning to return the modulus of the transforms and
    \code{TRUE} meaning to return the complex values themselves}
  \item{n}{second argument to \code{lpham()}}
  \item{...}{Ignored}
}
\details{
  \itemize{
    \item Function \code{is.near_integer(i)} returns \code{TRUE} if
    \code{i} is \dQuote{near} [that is, within \code{tol}] an integer;
    if the option is unset then \code{1e-11} is used.
    \item Function \code{is.nonpos()} returns \code{TRUE} if \code{i} is
    near a nonpositive integer
    \item Function \code{is.zero()} returns \code{TRUE} if \code{i} is,
    er, near zero
    \item Function \code{isgood()} checks for all elements of \code{x}
    having absolute values less than \code{tol}
    \item Function \code{thingfun()} transforms input vector \code{z} by
    each of the six members of the anharmonic group, viewed as a
    subgroup of the Mobius group of functions.  It returns a real
    six-column matrix with columns being the modulus of
    \eqn{z,z/(z-1),1-z,1/z,1/(1-z),1-1/z}.  These six columns
    correspond to the primary argument in equations 15.3.3 to 15.3.9,
    p551 of AMS-55
    \item Function \code{crit()} returns the two critical points,
    \eqn{\frac{1}{2}\pm\frac{\sqrt{3}i}{2}}{1/2 +/- sqrt(3)i/2}.  These
    points have unit modulus as do their six transforms by
    \code{thingfun()}
    \item Function \code{lpham()} returns the log of the Pochhammer
    function
    \eqn{log\left(\Gamma(x+n)/\Gamma(x)\right)}{log(Gamma(x+n)/Gamma(x))}
  }
}
  \author{Robin K. S. Hankin}
\note{
  Function \code{isgood()} uses zero as the default tolerance (argument
  \code{tol} passed in from \code{hypergeo()});
  compare the different meaning of \code{tol} used in
  \code{is.near_integer()}.

  Here, \dQuote{integer} means one of the sequence \eqn{0,\pm 1,\pm
  2,\ldots}{0, +/-1, +/-2, ...} [ie \emph{not} the Gaussian integers].
  
}
\examples{

is.near_integer(-3)

is.zero(4)
}
\keyword{math}

