\name{jacobian}
\alias{jacobian}
\alias{jacobian.default}
\title{Gradient of a Vector Valued Function}
\description{
Calculate the m by n numerical approximation of the gradient of a real
m-vector valued function with n-vector argument.
}
\usage{
    jacobian(func, x, method="Richardson", side=NULL, method.args=list(), ...) 

    \method{jacobian}{default}(func, x, method="Richardson", side=NULL,
       method.args=list(), ...)
}
\arguments{
    \item{func}{a function with a real (vector) result.}
    \item{x}{a real or real vector argument to func, indicating the point 
       at which the gradient is to be calculated.}
    \item{method}{one of \code{"Richardson"}, \code{"simple"}, or 
       \code{"complex"} indicating the method to use for the approximation.}
    \item{method.args}{arguments passed to method. See \code{\link{grad}}. 
       (Arguments not specified remain with their default values.)}
    \item{...}{any additional arguments passed to \code{func}.
          WARNING: None of these should have names matching other arguments of this function.}
    \item{side}{an indication of whether one-sided derivatives should be
       attempted (see details in function \code{\link{grad}}).}
}
\value{A real m by n matrix.}
\details{
   For \eqn{f:R^n -> R^m}{f:R^n -> R^m} calculate the \eqn{m x n}{m x n} 
   Jacobian \eqn{dy/dx}{dy/dx}.
   The function \code{jacobian} calculates a numerical approximation of the 
   first derivative of \code{func} at the point \code{x}. Any additional 
   arguments in \dots are also passed to \code{func}, but the gradient is not
   calculated with respect to these additional arguments.

   If method is "Richardson", the calculation is done by 
   Richardson's extrapolation. See \code{link{grad}} for more details.
   For this method  \code{method.args=list(eps=1e-4, d=0.0001, 
   zero.tol=sqrt(.Machine$double.eps/7e-7), r=4, v=2, show.details=FALSE)} 
   is set as the default.

   If method is "simple", the calculation is done using a simple epsilon
   difference. 
   For method "simple" \code{method.args=list(eps=1e-4)} is the
   default. Only \code{eps} is used by this method.

   If method is "complex", the calculation is done using the complex step
   derivative approach. See addition comments in \code{\link{grad}} before
   choosing this method.  
   For method "complex", \code{method.args} is ignored.
   The algorithm uses an \code{eps} of \code{.Machine$double.eps} which cannot
   (and should not) be modified. 
}
\seealso{
  \code{\link{grad}},
  \code{\link{hessian}},
  \code{\link[stats]{numericDeriv}}
  }
\examples{
   func2 <- function(x) c(sin(x), cos(x))
   x <- (0:1)*2*pi
   jacobian(func2, x)
   jacobian(func2, x, "complex")
}
\keyword{multivariate}
