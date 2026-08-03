\name{numericGradient}
\alias{numericGradient}
\alias{numericHessian}
\alias{numericNHessian}

\title{Functions to Calculate Numeric Derivatives}
\description{
  Calculate (central) numeric gradient and Hessian, including of
  vector-valued functions.
}
\usage{
numericGradient(f, t0, eps=1e-06, fixed, \dots)
numericHessian(f, grad=NULL, t0, eps=1e-06, fixed, \dots)
numericNHessian(f, t0, eps=1e-6, fixed, \dots)

}

\arguments{
  \item{f}{function to be differentiated.  The first argument must be
    the parameter vector with respect to which it is differentiated.
  For numeric gradient, \code{f} may return a (numeric) vector, for Hessian it
  should return a numeric scalar}
  \item{grad}{function, gradient of \code{f}}
  \item{t0}{vector, the parameter values}
  \item{eps}{numeric, the step for numeric differentiation}
  \item{fixed}{logical index vector, fixed parameters.
    Derivative is calculated only with respect to the parameters
    for which \code{fixed == FALSE}, \code{NA} is returned for the fixed
    parameters.  If
    missing, all parameters are treated as active.}
  \item{\dots}{furter arguments for \code{f}}
}
\details{
  \code{numericGradient} numerically differentiates a (vector valued)
  function with respect to it's (vector valued) argument.  If the
  functions value is a \eqn{N_{val} \times 1}{\code{N_val * 1}}
  vector and the argument is
  \eqn{N_{par} \times 1}{\code{N_par * 1}} vector, the resulting
  gradient
  is a \eqn{N_{val} \times N_{par}}{\code{NVal * NPar}}
  matrix. 

  \code{numericHessian} checks whether a gradient function is present.
  If yes, it calculates the gradient of the gradient, if not, it
  calculates the full
  numeric Hessian (\code{numericNHessian}).
}
\value{
  Matrix.  For \code{numericGradient}, the number of rows is equal to the
  length of the function value vector, and the number of columns is
  equal to the length of the parameter vector.

  For the \code{numericHessian}, both numer of rows and columns is
  equal to the length of the parameter vector.
}
\section{Warning}{
  Be careful when using numerical differentiation in optimization
  routines.  Although quite precise in simple cases, they may work very
  poorly in more complicated conditions.
  }
\author{Ott Toomet}
\seealso{\code{\link{compareDerivatives}}, \code{\link{deriv}}}

\examples{
# A simple example with Gaussian bell surface
f0 <- function(t0) exp(-t0[1]^2 - t0[2]^2)
numericGradient(f0, c(1,2))
numericHessian(f0, t0=c(1,2))

# An example with the analytic gradient
gradf0 <- function(t0) -2*t0*f0(t0)
numericHessian(f0, gradf0, t0=c(1,2))
# The results should be similar as in the previous case

# The central numeric derivatives are often quite precise
compareDerivatives(f0, gradf0, t0=1:2)
# The difference is around 1e-10
}
\keyword{math}
\keyword{utilities}
