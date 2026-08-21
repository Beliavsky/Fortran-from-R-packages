\name{grad}
\alias{grad}
\alias{grad.default}
\title{Numerical Gradient of a Function}
\description{Calculate the gradient of a function by numerical approximation.}
\usage{
     grad(func, x, method="Richardson", side=NULL, method.args=list(), ...) 

    \method{grad}{default}(func, x, method="Richardson", side=NULL,
      method.args=list(), ...)
 }
\arguments{
    \item{func}{a function with a scalar real result (see details).}
    \item{x}{a real scalar or vector argument to func, indicating the 
      point(s) at which the gradient is to be calculated.}
    \item{method}{one of \code{"Richardson"}, \code{"simple"}, or 
       \code{"complex"} indicating the method to use for the approximation.}
    \item{method.args}{arguments passed to method. Arguments not specified 
       remain with their default values as specified in details}
    \item{side}{an indication of whether one-sided derivatives should be
       attempted (see details).}
    \item{...}{an additional arguments passed to \code{func}.
          WARNING: None of these should have names matching other arguments of this function.}
}
\value{A real scalar or vector of the approximated gradient(s).}

\details{
   The function \code{grad} calculates a numerical approximation of the 
   first derivative of \code{func} at the point \code{x}. Any additional 
   arguments in \dots are also passed to \code{func}, but the gradient is not
   calculated with respect to these additional arguments.
   It is assumed \code{func}  is a scalar value function. If a vector \code{x} 
   produces a scalar
   result then \code{grad} returns the numerical approximation of the gradient
   at the point \code{x} (which has the same length as \code{x}).
   If a vector \code{x} produces a vector result then the result must have the
   same length as \code{x}, and it is assumed that this corresponds to applying
   the function to each of its arguments (for example, \code{sin(x)}). 
   In this case \code{grad} returns the
   gradient at each of the points in \code{x} (which also has the same length 
   as \code{x} -- so be careful). An alternative for vector valued functions is
   provided by \code{\link{jacobian}}.
  
   If method is "simple", the calculation is done using a simple epsilon
   difference. 
   For method "simple" \code{method.args=list(eps=1e-4)} is the
   default. Only \code{eps} is used by this method.

   If method is "complex", the calculation is done using the complex step
   derivative approach of Lyness and Moler, described in  Squire and Trapp. 
   This method requires that the function be able to handle complex valued 
   arguments and return the appropriate complex valued result, 
   even though the user may only be interested in the real-valued derivatives. 
   It also requires that the complex function be analytic. (This might be thought 
   of as the complex equivalent of the requirement for continuity and smoothness 
   of a real valued function.) 
   So, while this method is extremely powerful it is applicable to
   a very restricted class of functions. \emph{Avoid this method if you do not 
   know that your function is suitable. Your mistake may not be caught and the
   results will be spurious.}
   For cases where it can be used,
   it is faster than Richardson's extrapolation, and
   it also provides gradients that are correct to machine precision (16 digits).   
   For method "complex", \code{method.args} is ignored.
   The algorithm uses an \code{eps} of \code{.Machine$double.eps} which cannot
   (and should not) be modified. 

   If method is "Richardson", the calculation
   is done by Richardson's extrapolation (see e.g. Linfield and Penny, 1989,
   or Fornberg and Sloan, 1994.)
   This method should be used if accuracy, as opposed to speed, is important
   (but see method "complex" above). 
   For this method 
   \code{method.args=list(eps=1e-4, d=0.0001, zero.tol=sqrt(.Machine$double.eps/7e-7), 
   r=4, v=2, show.details=FALSE)} is set as the default.
    \code{d} gives the fraction of \code{x} to use for the initial numerical 
      approximation. The default means the initial approximation uses
      \code{0.0001 * x}.
    \code{eps} is used instead of \code{d} for elements of \code{x} which are 
       zero (absolute value less than zero.tol).
    \code{zero.tol} tolerance used for deciding which elements of \code{x} are 
       zero.
    \code{r} gives the number of Richardson improvement iterations (repetitions
       with successly smaller \code{d}. The default \code{4} general provides 
       good results, but this can be increased to \code{6} for improved
       accuracy at the cost of more evaluations.
    \code{v} gives the reduction factor.
    \code{show.details} is a logical indicating if detailed calculations should 
        be shown.
   
  The general approach in the Richardson method is to iterate for \code{r} 
  iterations from initial 
  values for interval value \code{d},  using reduced factor \code{v}.
  The the first order approximation to the derivative with respect 
  to \eqn{x_{i}}{x_{i}} is

      \deqn{f'_{i}(x) = <f(x_{1},\dots,x_{i}+d,\dots,x_{n}) -
               f(x_{1},\dots,x_{i}-d,\dots,x_{n})>/(2*d)}{%
	    f'_{i}(x) = <f(x_{1},\dots,x_{i}+d,\dots,x_{n}) -
               f(x_{1},\dots,x_{i}-d,\dots,x_{n})>/(2*d)}
       
  This is repeated \code{r} times  with successively smaller \code{d}  and 
  then Richardson extraplolation is applied.
  
  If elements of \code{x} are near zero the multiplicative interval calculation
  using \code{d} does not work, and for these elements an additive calculation
  using \code{eps} is done instead. The argument \code{zero.tol} is used
  determine if an element should be considered too close to zero. 
  In the iterations, interval is successively reduced to eventual 
  be \code{d/v^r} and the square of this value is used in second derivative 
  calculations (see \code{\link{genD}}) so the 
  default \code{zero.tol=sqrt(.Machine$double.eps/7e-7)} is set to ensure the
  interval is bigger than \code{.Machine$double.eps} with the default \code{d},
  \code{r}, and \code{v}.
  
  If \code{side} is \code{NULL} then it is assumed that the point at which the
  calculation is being done is interior to the domain of the function. If the
  point is on the boundary of the domain then \code{side} can be used to 
  indicate which side of the point \code{x} should be used for the calculation.
  If not \code{NULL} then it should be a vector of the same length as \code{x}
  and have values \code{NA}, \code{+1}, or \code{-1}. \code{NA} indicates that
  the usual calculation will be done, while \code{+1}, or \code{-1} indicate
  adding or subtracting from the parameter point \code{x}. The argument
  \code{side} is not supported for all methods.
  
  Since usual calculation with method "simple" uses only a small \code{eps} 
  step to one side, the only effect of argument \code{side} is to determine the
  direction of the step. The usual calculation with method "Richardson" is 
  symmetric, using steps to both sides. The effect of argument \code{side} 
  is to take a double sized step to one side, and no step to the other side.
  This means that the center of the Richardson extrapolation steps is moving
  slightly in the reduction, and is not exactly on the boundary. 
  (Warning: I am not aware of theory or published
  experimental evidence to support this, but the results in my limited testing
  seem good.)
}
\references{ 
   Linfield, G. R. and Penny, J. E. T. (1989) \emph{Microcomputers in Numerical 
   Analysis}. New York: Halsted Press.

  Fornberg, B. and Sloan, D, M. (1994) ``A review of pseudospectral methods 
  for solving partial differential equations.'' \emph{Acta Numerica}, 3, 203-267.

  Lyness, J. N. and Moler, C. B. (1967) ``Numerical Differentiation of Analytic 
   Functions.'' \emph{SIAM Journal for Numerical Analysis},
    4(2), 202-210.

  Squire, William and Trapp, George (1998) ``Using Complex Variables to Estimate
   Derivatives of Real Functions.'' \emph{SIAM Rev},
    40(1), 110-112.
}
\seealso{
  \code{\link{jacobian}},
  \code{\link{hessian}},
  \code{\link{genD}},
  \code{\link[stats]{numericDeriv}}
   }
\examples{
  grad(sin, pi)
  grad(sin, (0:10)*2*pi/10)
  func0 <- function(x){ sum(sin(x))  }
  grad(func0 , (0:10)*2*pi/10)

  func1 <- function(x){ sin(10*x) - exp(-x) }

  curve(func1,from=0,to=5)

  x <- 2.04
  numd1 <- grad(func1, x)
  exact <- 10*cos(10*x) + exp(-x)
  c(numd1, exact, (numd1 - exact)/exact)

  x <- c(1:10)
  numd1 <- grad(func1, x)
  numd2 <- grad(func1, x, "complex")
  exact <- 10*cos(10*x) + exp(-x)
  cbind(numd1, numd2, exact, (numd1 - exact)/exact, (numd2 - exact)/exact)

  sc2.f <- function(x){
    n <- length(x)
    sum((1:n) * (exp(x) - x)) / n
    }

  sc2.g <- function(x){
    n <- length(x)
    (1:n) * (exp(x) - 1) / n
    }

  x0 <- rnorm(100)
  exact <- sc2.g(x0)

  g <- grad(func=sc2.f, x=x0)
  max(abs(exact - g)/(1 + abs(exact)))

  gc <- grad(func=sc2.f, x=x0, method="complex")
  max(abs(exact - gc)/(1 + abs(exact)))

  f <- function(x) if(x[1]<=0) sum(sin(x)) else  NA
  grad(f, x=c(0,0), method="Richardson", side=c(-1,  1))
}
\keyword{multivariate}
