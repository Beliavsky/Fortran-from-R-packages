\name{nnls-package}
\alias{nnls-package}
\docType{package}
\title{The Lawson-Hanson NNLS implementation of non-negative least squares}
\description{
  An R interface to the Lawson-Hanson
  NNLS implementation of an algorithm
  for non-negative linear least squares 
  that solves the least squares problem
  \eqn{\min{\parallel A x = b \parallel_2}}
  with the constraint \eqn{x \ge 0} where
  \eqn{x \in R^n, b \in R^m}  and \eqn{A} is an
  \eqn{m \times n} matrix. 
  Also allows the combination of non-negative and non-positive
  constraints on \eqn{x}. 
}

\references{
Lawson CL, Hanson RJ (1974). Solving Least Squares Problems. Prentice
Hall, Englewood Cliffs, NJ.

Lawson CL, Hanson RJ (1995). Solving Least Squares Problems. Classics
in Applied Mathematics. SIAM, Philadelphia.
}

\keyword{ package }
\seealso{ \link{nnls}, \link{nnnpls},
  the method \code{"L-BFGS-B"} for \link{optim},
   \link[quadprog]{solve.QP}, \link[bvls]{bvls}
} 

