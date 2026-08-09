\name{nlsDeriv}
\alias{nlsDeriv}
\alias{codeDeriv}
\alias{fnDeriv}
\title{nlsDeriv
  Functions to take symbolic derivatives.}
\usage{
nlsDeriv(expr, name, derivEnv = sysDerivs, do_substitute = FALSE, verbose = FALSE, ...)

codeDeriv(expr, namevec, hessian = FALSE, derivEnv = sysDerivs, 
          do_substitute = FALSE, verbose = FALSE, ...)

fnDeriv(expr, namevec, args = all.vars(expr), env = environment(expr), 
          do_substitute = FALSE, verbose = FALSE, ...)
}
\arguments{
\item{expr}{An expression represented in a variety of ways. See Details.}

\item{name}{The name of the variable with respect to which the derivative will be computed.}

\item{derivEnv}{The environment in which derivatives are stored.}

\item{do_substitute}{If \code{TRUE}, use \code{\link{substitute}} to get the expression passed as
\code{expr}, otherwise evaluate it.}

\item{verbose}{If \code{TRUE}, then diagnostic output will be printed as derivatives
and simplifications are recognized.}

\item{...}{Additional parameters which will be passed to \code{codeDeriv}
from \code{fnDeriv}, and to \code{nlsSimplify} from 
\code{nlsDeriv} and \code{codeDeriv}.}

\item{namevec}{Character vector giving the variable names with respect to 
which the derivatives will be taken.}

\item{hessian}{Logical indicator of whether the 2nd derivatives should also be computed.}

\item{args}{Desired arguments for the function.  See Details below.}

\item{env}{The environment to be attached to the created function.
If \code{NULL}, the caller's frame is used.}
}
\description{
Compute derivatives of simple expressions symbolically, allowing user-specified derivatives.
}
\details{   
 Functions \code{nlsDeriv} and \code{codeDeriv} are designed as replacements 
 for the \pkg{stats} package functions \code{\link{D}} and \code{\link{deriv}}
 respectively, though the argument lists do not match exactly.
 
 The \code{nlsDeriv} function computes a symbolic derivative of an expression
 or language object.  Known derivatives are stored in
 \code{derivEnv}; the default \code{sysDerivs} contains expressions for
 all of the derivatives recognized by \code{\link{deriv}}, but in
 addition allows differentiation with respect to any parameter
 where it makes sense.  It also allows the derivative of \code{abs}
 and \code{sign}, using an arbitrary choice of 0 at the discontinuities.
 
 The \code{codeDeriv} function computes
 an expression for efficient calculation of the expression value together
 with its gradient and optionally the Hessian matrix.
 
 The \code{fnDeriv} function wraps the \code{codeDeriv} result
 in a function.  If the \code{args} are given as a character
 vector (the default), the arguments will have those names,
 with no default values.  Alternatively, a custom argument list with default values can
 be created using \code{\link{alist}}; see the example below.
 
 The \code{expr} argument will be converted to a
 language object using \code{\link{dex}} (but note
 the different default for \code{do_substitute}).  
 Normally it should be a formula with no left
 hand side, e.g. \code{ ~ x^2 }, or an expression vector
 e.g. \code{ expression(x, x^2, x^3) }, or a language
 object e.g. \code{quote(x^2)}.  In \code{codeDeriv} and
 \code{fnDeriv} the expression vector must be of length 1.
 
 The \code{newDeriv} function is used to define a new derivative.
 The \code{expr} argument should match the header of the function as a
 call to it (e.g. as in the help pages), and the \code{deriv} argument
 should be an expression giving the derivative, including calls to
 \code{D(arg)}, which will not be evaluated, but will be substituted
 with partial derivatives of that argument with respect to \code{name}.
 See the examples below.  
 
 If \code{expr} or \code{deriv} is missing in a call to
 \code{newDeriv()}, it will return the currently saved derivative
 record from \code{derivEnv}.  If \code{name} is missing in a call to
 \code{nlsDeriv} with a function call, it will print a message describing
 the derivative formula and return \code{NULL}.
 
 To handle functions which act differently if a parameter is
 missing, code the default value of that parameter to \code{.MissingVal},
 and give a derivative that is conditional on \code{missing()}
 applied to that parameter.  See the derivatives of \code{"-"} and \code{"+"} 
 in the file \code{derivs.R} for an example.
}
 
\value{
 If \code{expr} is an expression vector, \code{nlsDeriv} and \code{nlsSimplify}
 return expression vectors containing the response.  
 For formulas or language objects, a language object is returned.
 
 \code{codeDeriv} always returns a language object.
 
 \code{fnDeriv} returns a closure (i.e. a function).
 
 \code{nlsDeriv} returns the symbolic derivative of the expression.
 
 \code{newDeriv} with \code{expr} and \code{deriv} specified is
 called for the side effect of recording the derivative in \code{derivEnv}.
 If \code{expr} is missing, it will return the list of names of functions
 for which derivatives are recorded.  If \code{deriv} is missing, it
 will return its record for the specified function.
}
\note{
 \code{newDeriv(expr, deriv, ...)} will issue a warning
 if a different definition for the derivative exists
 in the derivative table.
}
\author{
  Duncan Murdoch
}
\seealso{\code{\link{deriv}}}
 
\examples{
 nlsDeriv(~ sin(x+y), "x")

 f <- function(x) x^2
 newDeriv(f(x), 2*x*D(x))
 nlsDeriv(~ f(abs(x)), "x")
 
 nlsDeriv(~ pnorm(x, sd=2, log = TRUE), "x")
 fnDeriv(~ pnorm(x, sd = sd, log = TRUE), "x")
 f <- fnDeriv(~ pnorm(x, sd = sd, log = TRUE), "x", args = alist(x =, sd = 2))
 f
 f(1)
 100*(f(1.01) - f(1))  # Should be close to the gradient
 
       # The attached gradient attribute (from f(1.01)) is
       # meaningless after the subtraction.
       
 # Multiple point example
 xvals <- c(1, 3, 4.123)
 print(f(xvals))
 # Getting a hessian matrix
 f2 <- ~ (x-2)^3*y - y^2
 mydf2 <- fnDeriv(f2, c("x","y"), hessian=TRUE)
 # display the resulting function
 print(mydf2)
 x <- c(1, 2)
 y <- c(0.5, 0.1)
 evalmydf2 <- mydf2(x, y)
 print(evalmydf2)
 # the first index of the hessian attribute is the point at which we want the hessian
 hmat1 <- as.matrix(attr(evalmydf2,"hessian")[1,,])
 print(hmat1)
 hmat2 <- as.matrix(attr(evalmydf2,"hessian")[2,,])
 print(hmat2)
}
\keyword{
   jacobian
}
