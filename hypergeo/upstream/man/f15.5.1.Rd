\name{f15.5.1}
\alias{f15.5.1}
\alias{f15.5.1}
\alias{hypergeo_press}
\alias{hypergeo_func}
\alias{to_real}
\alias{to_complex}
\alias{complex_ode}
\alias{semicircle}
\alias{semidash}
\alias{straight}
\alias{straightdash}
\title{
Hypergeometric functions via direct numerical integration
}
\description{
  The hypergeometric function may be evaluated using Gauss's
  differential equation 15.5.1: 
  \deqn{z(1-z)w''+(c-(a+b+1)z)w'-abw=0}

using a start value away from the three singular points.  This page
documents a suite of related functionality.

}
\usage{
hypergeo_press(A,B,C,z, ...)
f15.5.1(A, B, C, z, startz, u, udash, give=FALSE, ...)
hypergeo_func(Time, State, Pars, u, udash)
to_real(o)
to_complex(p)
complex_ode(y, times, func, parms=NA, method=NULL, u, udash, ...)
semicircle(t,z0,z1,clockwise=TRUE)
semidash(t,z0,z1,clockwise=TRUE)
straight(t,z0,z1)
straightdash(t,z0,z1)
}
\arguments{
  \item{A,B,C,z}{Standard parameters for the hypergeometric function}
  \item{u,udash}{Functions to specify the path of integration, and its
    derivative}
  \item{give}{In function \code{f15.5.1()}, Boolean with \code{TRUE}
    meaning to return extra information from \code{ode()} and default
    \code{FALSE} meaning to return only the evaluated function}
  \item{startz}{In function \code{f15.5.1(}, the start position of the path}
  \item{\dots}{Further arguments passed to \code{ode()}}
  \item{o,p}{Real and complex objects to be coerced to each other in
    \code{to_real()} and \code{to_complex()}}
  \item{y, times, func, parms, method}{In function \code{complex_ode()},
    arguments matching those of \code{ode()}}
  \item{t,z0,z1,clockwise}{Arguments for the standard path functions
      \code{semicircle()} et seq: \code{u} is the primary argument
      (real, \eqn{0\leq u\leq 1}{1<=u<=1}); \code{z0} and \code{z1} are
      the start and end points of the path; and \code{clockwise} is
      Boolean, indicating whether the path proceeds clockwise or not}
  \item{Time, State, Pars}{arguments matchin those of standard examples
    in the \code{deSolve} package}
}
\details{

  Function \code{hypergeo_press()} is the most user-friendly of the
  functions documented here.  It performs integration of Gauss's ODE,
  along a straight line path from the start-point to \code{z}.  It follows
  Press et al's suggestion of start-point.
  
  Function \code{f15.5.1()} is a little more flexible in that it allows
  the user to choose a start point and an integration path.

  Function \code{complex_ode()} is a complex generalization of
  \code{ode()} of package \pkg{deSolve}; function \code{hypergeo_func}
  is an internal function, designed for use with \code{complex_ode()},
  that specifies the Gauss ODE which is satisified by the hypergeometric
  function.

  Functions \code{to_real()} and \code{to_complex()} are internal
  functions which coerce from real to complex and back; they are needed
  because \code{ode()} deals only with real values.

  Functions \code{semicircle()} and \code{straight()} are helper
  functions which specify straight or semicircular paths from \code{z0}
  to \code{z1}; note that \code{f(0)=z0} and \code{f(1)=z1}.  Functions
  \code{semidash()} and \code{straightdash()} provide the differential
  of the path.
}
\references{
W. H. Press et al.  1997. \emph{Numerical Recipes in C}.  Cambridge
  University Press, Second Edition.
}
\author{
Robin K. S. Hankin
}
\note{
Accuracy is low compared with the other methods in the package.
}

\seealso{
\code{\link{hypergeo_residue}}
}
\examples{

hypergeo_press(A=pi,B=sqrt(2),C=1.4,z=1-2i)
hypergeo      (A=pi,B=sqrt(2),C=1.4,z=1-2i)


jj1 <- 
f15.5.1(
    A=1.1, B=2.2, C=3.3, z=3+0.5i, startz=0.5,
        u    =function(u){semicircle(u,0.5,3+0.5i,FALSE)},
        udash=function(u){semidash(u,0.5,3+0.5i,FALSE)}
        )

jj2 <-
f15.5.1(
    A=1.1, B=2.2, C=3.3, z=3+0.5i, startz=0.5,
        u    =function(u){semicircle(u,0.5,3+0.5i,TRUE)},
        udash=function(u){semidash(u,0.5,3+0.5i,TRUE)}
        )



jj3 <- hypergeo(    A=1.1, B=2.2, C=3.3, z=3+0.5i)
## First one agrees with jj3=hypergeo(...), the second one does not 


# Now try the Airy Ai function;  satisfies f'' =  zf:

pars <- c(kay = 1+1i, ell = 0.1+0.2i)  # not actually used
airy_ai_func <- function(Time, State, Pars, u, udash) {
    with(as.list(c(to_complex(State), to_complex(Pars))), {

      z <- u(Time)
      dz <- udash(Time)
        
      dF <- Fdash*dz
      dFdash <-  z*F*dz # could use kay and ell from pars here if necessary
        
        ## coerce back to real:
        out <- to_real(c(dF,dFdash))
        names(out) <- names(State)
        return(list(out))
    })
}

complex_ode(
    y     = c(F = 1/3^(2/3)/gamma(2/3), Fdash= -1/3^(1/3)/gamma(1/3)),
    times = seq(0,1,by=0.1),
    func  = airy_ai_func,
    parms = pars,
    u     = function(t){straight(t,0,1)},
    udash = function(t){straightdash(t,0,1)}
)

# Look at the last line for the value at 1.
# compare gsl: Ai(1) = 0.1352924 ; Ai'(1) = -0.1591474

# ...although in this case there is actually a hypergeometric series
#  representation:

f <- function(z){
    return(
        +genhypergeo(U=NULL,L=2/3,z^3/9)/3^(2/3)/gamma(2/3)
        -genhypergeo(U=NULL,L=4/3,z^3/9)/3^(1/3)/gamma(1/3)*z
        ) 
}

f(1)

}
