\name{limSolve-package}
\alias{limSolve-package}
\alias{limSolve}
\docType{package}
\title{
  Solving Linear Inverse Models
}
\description{
  Functions that:

  (1.) Find the minimum/maximum of a linear or quadratic function:
   min or max (f(x)), where \eqn{f(x) = ||Ax-b||^2} or \eqn{f(x) = sum(ai*xi)}
   subject to equality constraints \eqn{Ex=f} and/or inequality constraints
   \eqn{Gx>=h}.

  (2.) Sample an underdetermined- or overdetermined system \eqn{Ex=f}
  subject to \eqn{Gx>=h}, and if applicable \eqn{Ax~=b}.

  (3.) Solve a linear system \eqn{Ax=B} for the unknown x.
  Includes banded and tridiagonal linear systems.

  The package calls Fortran functions from LINPACK
}

\references{
Van den Meersche K, Soetaert K, Van Oevelen D (2009). xsample(): An R Function for
Sampling Linear Inverse Problems.
Journal of Statistical Software, Code Snippets, 30(1), 1-15.

}

\details{

  limSolve is designed for solving linear inverse models (LIM).

  These consist of linear equality and, or inequality conditions,
  which can be solved either by least squares or by linear programming
  techniques.

  Amongst the possible applications are: food web quantification,
  flux balance analysis (e.g. metabolic networks),
  compositional estimation, and operations research problems.

  The package contains several examples to exemplify its use
}

\author{
  Karline Soetaert (Maintainer),
  
  Karel Van den Meersche

  Dick van Oevelen
}
%\references{
%% to be filled
%}
\seealso{
 \code{\link{Blending}}, \code{\link{Chemtax}}, \code{\link{RigaWeb}},
 \code{\link{E_coli}}, \code{\link{Minkdiet}} the examples.

 \code{\link{ldei}}, \code{\link{lsei}},\code{\link{linp}}, \code{\link{ldp}},
 \code{\link{nnls}} to solve LIM

 \code{\link{xranges}}, \code{\link{varranges}}  to estimate ranges of
 unknowns and variables

 \code{\link{xsample}}, \code{\link{varsample}} to create a random sample
 of unknowns and variables

 \code{\link{Solve}}, \code{\link{Solve.banded}}, \code{\link{Solve.tridiag}},
 to solve non-square, banded and tridiagonal linear systems of equations.
 
 \code{\link{resolution}} row and column resolution of a matrix
 
 package vignette \code{limSolve}
}

\examples{
\dontrun{
## show examples (see respective help pages for details)
example(Blending)
example(Chemtax)
example(E_coli)
example(Minkdiet)

## run demos
demo("limSolve")

## open the directory with original E_coli input file
browseURL(paste(system.file(package="limSolve"), "/doc", sep=""))

## show package vignette with tutorial about xsample
vignette("xsample")

## show main package vignette
vignette("limSolve")
}
}

\keyword{ package }