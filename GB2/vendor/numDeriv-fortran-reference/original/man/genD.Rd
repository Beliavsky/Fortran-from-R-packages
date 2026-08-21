\name{genD}
\alias{genD}
\alias{genD.default}

\title{Generate Bates and Watts D Matrix}
\description{Generate a matrix of function derivative information.}
\usage{
    genD(func, x, method="Richardson",
                   method.args=list(), ...)
    \method{genD}{default}(func, x, method="Richardson",
      method.args=list(), ...)
}
\arguments{
    \item{func}{a function for which the first (vector) argument 
        is used as a parameter vector.}
    \item{x}{The parameter vector first argument to \code{func}.}
    \item{method}{one of \code{"Richardson"} or \code{"simple"} indicating 
    the method to use for the aproximation.}
    \item{method.args}{arguments passed to method.  See \code{\link{grad}}. 
       (Arguments not specified remain with their default values.)}
    \item{...}{any additional arguments passed to \code{func}.
          WARNING: None of these should have names matching other arguments of this function.}
}
\value{
A list with elements as follows:
   \code{D} is a matrix of first and second order partial
      derivatives organized in the same manner as Bates and 
      Watts, the number of rows is equal to the length of the result of
      \code{func}, the first p columns are the Jacobian, and the 
      next p(p+1)/2 columns are the lower triangle of the second derivative
      (which is the Hessian for a scalar valued \code{func}).
   \code{p} is the length of \code{x} (dimension of the parameter space).
   \code{f0} is the function value at the point where the matrix \code{D} 
        was calculated. 
   The  \code{genD} arguments \code{func}, \code{x}, \code{d}, \code{method},
   and  \code{method.args} also are returned in the list.
}
\details{
   The derivatives are calculated numerically using Richardson improvement.
   Methods "simple" and "complex" are not supported in this function.
   The "Richardson" method calculates a numerical approximation of the first 
   and second derivatives of \code{func} at the point \code{x}. 
   For a scalar valued function these are the gradient vector and 
   Hessian matrix. (See \code{\link{grad}} and \code{\link{hessian}}.)
   For a vector valued function the first derivative is the Jacobian matrix 
   (see \code{\link{jacobian}}). 
   For the Richardson method 
   \code{method.args=list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), 
   r=4, v=2)} is set as the default.
   See \code{\link{grad}}
   for more details on the Richardson's extrapolation parameters. 

   A simple approximation to the first order derivative with respect 
   to \eqn{x_i}{x_i} is 

   \deqn{f'_{i}(x) = <f(x_{1},\dots,x_{i}+d,\dots,x_{n}) -
               f(x_{1},\dots,x_{i}-d,\dots,x_{n})>/(2*d)}{%
 	 f'_{i}(x) = <f(x_{1},\dots,x_{i}+d,\dots,x_{n}) -
               f(x_{1},\dots,x_{i}-d,\dots,x_{n})>/(2*d)}

  A simple approximation to the second order derivative with respect 
  to \eqn{x_i}{x_i} is 

  \deqn{f''_{i}(x) = <f(x_{1},\dots,x_{i}+d,\dots,x_{n}) -
                   2 *f(x_{1},\dots,x_{n}) +
                    f(x_{1},\dots,x_{i}-d,\dots,x_{n})>/(d^2) }{%
	f''_{i}(x) = <f(x_{1},\dots,x_{i}+d,\dots,x_{n}) -
                   2 *f(x_{1},\dots,x_{n}) +
                    f(x_{1},\dots,x_{i}-d,\dots,x_{n})>/(d^2) }	    

   The second order derivative with respect to \eqn{x_i, x_j}{x_i, x_j} is 

  \deqn{f''_{i,j}(x) = <f(x_{1},\dots,x_{i}+d,\dots,x_{j}+d,\dots,x_{n}) -
                    2 *f(x_{1},\dots,x_{n}) + }{%
        f''_{i,j}(x) = <f(x_{1},\dots,x_{i}+d,\dots,x_{j}+d,\dots,x_{n}) -
                    2 *f(x_{1},\dots,x_{n}) + }

  \deqn{f(x_{1},\dots,x_{i}-d,\dots,x_{j}-d,\dots,x_{n})>/(2*d^2) -
		      (f''_{i}(x) + f''_{j}(x))/2 }{%
        f(x_{1},\dots,x_{i}-d,\dots,x_{j}-d,\dots,x_{n})>/(2*d^2) -
		      (f''_{i}(x) + f''_{j}(x))/2 }


   Richardson's extrapolation is based on these formula with the \code{d} 
   being reduced in the extrapolation iterations. In the code, \code{d} is
   scaled to accommodate parameters of different magnitudes. 

   \code{genD} does \code{1 + r (N^2 + N)} evaluations of the function
   \code{f}, where \code{N} is the length of \code{x}.
}

\references{ 
   Linfield, G.R. and Penny, J.E.T. (1989) "Microcomputers in Numerical Analysis."
    Halsted Press.
   
   Bates, D.M. & Watts, D. (1980), "Relative Curvature Measures of Nonlinearity."
   J. Royal Statistics Soc. series B, 42:1-25

   Bates, D.M. and Watts, D. (1988) "Non-linear Regression Analysis and Its Applications."
   Wiley.
}

\seealso{
   \code{\link{hessian}}, 
   \code{\link{grad}}
}
\examples{
    func <- function(x){c(x[1], x[1], x[2]^2)}
    z <- genD(func, c(2,2,5))
}
\keyword{multivariate}
