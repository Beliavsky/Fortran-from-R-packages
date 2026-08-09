\name{nlsr.package}
\alias{nlsr.package}
\title{nlsr-package
 
Tools for solving nonlinear least squares problems
 
The package provides some tools related to using the Nash variant
of Marquardt's algorithm for nonlinear least squares. Jacobians can
usually be developed by automatic or symbolic derivatives.}
\usage{
nlsr.package()
}
\description{
nlsr-package
 
Tools for solving nonlinear least squares problems
 
The package provides some tools related to using the Nash variant
of Marquardt's algorithm for nonlinear least squares. Jacobians can
usually be developed by automatic or symbolic derivatives.
}
\details{
This package includes methods for solving nonlinear least squares problems
specified by a modeling expression and given a starting vector of named 
paramters. Note: You must provide an expression of the form
    lhs ~ rhsexpression
so that the residual expression 
    rhsexpression - lhs
can be computed. The expression can be enclosed in quotes, and this seems to give
fewer difficulties with R. Data variables must already be defined, either within 
the parent environment or else in the dot-arguments. Other symbolic elements in
the modeling expression must be standard functions or else parameters that are 
named in the start vector.

The main functions in \code{nlsr} are:

nlfb   Nash variant of the Marquardt procedure for nonlinear least squares,
   with bounds constraints, using a residual and optionally Jacobian
   described as \code{R} functions.
   
nlxb  Nash variant of the Marquardt procedure for nonlinear least squares,
  	with bounds constraints, using an expression to describe the residual via
  	an \code{R} modeling expression. The Jacobian is computed via symbolic
  	differentiation.
  	
wrapnlsr  Uses \code{nlxb} to solve nonlinear least squares then calls 
   \code{nls()} to create an object of type nls. \code{nlsr} is an alias
for \code{wrapnlsr}
   
model2rjfun  returns a function with header \code{function(prm)}, which
   evaluates the residuals (and if jacobian is TRUE the Jacobian matrix) 
   of the model at \code{prm}.  The residuals are defined to be the 
   right hand side of \code{modelformula} minus the left hand side.
   
model2ssgrfun  returns a function with header \code{function(prm)}, which
    evaluates the sum of squared residuals (and if gradient is \code{TRUE} the
    gradient vector) of the model at \code{prm}.
    
modelexpr  returns the expression used to calculate the vector of 
        residuals (and possibly the Jacobian) used in the previous functions.
}
\references{
Nash, J. C. (1979, 1990) _Compact Numerical Methods for Computers.
    Linear Algebra and Function Minimisation._ Adam Hilger./Institute
    of Physics Publications

  Nash, J. C. (2014) _Nonlinear Parameter Optimization Using R Tools._
    Wiley
}
\author{
John C Nash and Duncan Murdoch
}
\keyword{nls}
\keyword{nonlinear&least&squares}
