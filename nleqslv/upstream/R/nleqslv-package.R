# Copyright (c) 2026 Rob Carnell

#' The \pkg{nleqslv} package provides two algorithms for solving (dense)
#' nonlinear systems of equations.
#'
#' The methods provided are \itemize{ \item a
#' Broyden Secant method where the matrix of derivatives is updated after each
#' major iteration using the Broyden rank 1 update.  \item a full Newton method
#' where the Jacobian matrix of derivatives is recalculated at each iteration }
#' Both methods utilize global strategies such as line search or trust region
#' methods whenever the standard Newton/Broyden step does not lead to a point
#' closer to a root of the equation system. Both methods can also be used
#' without a norm reducing global strategy. Line search may be either cubic,
#' quadratic or geometric. The trust region methods are either the double
#' dogleg method, the Powell single dogleg method or a Levenberg-Marquardt type
#' method.
#'
#' There is a facility for specifying that the Jacobian is banded; this can
#' significantly speedup the calculation of a numerical Jacobian when the
#' number of sub- and super diagonals is small compared to the size of the
#' system of equations. For example the Jacobian of a tridiagonal system can be
#' calculated with only three evaluations of the function.
#'
#' The package provides an option to attempt to solve the system of equations
#' when the Jacobian is singular or ill-conditioned using an approximation to
#' the Moore-Penrose pseudoinverse of the Jacobian.
#'
#' The algorithms provided in this package are derived from Dennis and Schnabel
#' (1996). The code is written in Fortran 77 and Fortran 95 and uses Lapack and
#' BLAS routines as provided by the R system.
#'
#' @aliases nleqslv-package nleqslv.Intro
#' @references Dennis, J.E. Jr and Schnabel, R.B. (1996), \emph{Numerical
#' Methods for Unconstrained Optimization and Nonlinear Equations}, Siam.
#' @keywords internal
"_PACKAGE"
