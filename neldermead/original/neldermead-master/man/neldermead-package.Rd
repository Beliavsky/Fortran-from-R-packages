\name{neldermead-package}
\alias{neldermead-package}
\docType{package}
\title{
  R port of the Scilab neldermead module
}
\description{
  The goal of this package is to provide a Nelder-Mead direct search
  optimization method. That Nelder-Mead algorithm may be used in the following
  optimization context: \itemize{
    \item there is no need to provide the derivatives of the objective
      function,
    \item the number of parameters is small (up to 10-20),
    \item there are bounds and/or non linear constraints.
    }

  \strong{Design}

  This package provides the following components:\itemize{
    \item \pkg{neldermead} provides various Nelder-Mead variants and manages for
      Nelder-Mead specific settings, such as the method to compute the initial
      simplex, the specific termination criteria,
    \item \code{fminsearch} provides a simplified Nelder-Mead algorithm.
      Specific termination criteria, initial simplex and auxiliary settings are 
      automatically configured.
    \item \code{fminbnd} provides a simplified Box algorithm, ie the equivalent
      of \code{fminsearch} for unconstrained search.
    \item \code{optimset}, \code{optimget} provide commands to emulate their
      Scilab counterparts.
    \item \code{optimplotfunccount}, \code{optimplotx} and \code{optimplotfval}
      provide plotting features for the fminsearch function (Not implemented
      yet).
    \item \code{nmplot} provides a high-level component which provides directly
      output pictures for Nelder-Mead algorithm. (Not implemented yet).
    }
  The current component is based on the following packages \itemize{
    \item \pkg{optimbase}: provides an abstract class for a general optimization
      component, including the number of variables, the minimum and maximum
      bounds, the number of non linear inequality constraints, the loggin
      system, various termination criteria, the cost function, etc...
    \item \pkg{optimsimplex}: provides a class to manage a simplex made of an
      arbitrary number of vertices, including the computation of a simplex by
      various methods (axes, regular, Pfeffer's, randomized bounds), the
      computation of the size by various methods (diameter, sigma+, sigma-,
      etc...),
    }
  
  \strong{Features}

  The following is a list of features the Nelder-Mead prototype algorithm
  currently provides:\itemize{
    \item Provides 3 algorithms, including \itemize{
      \item the fixed shape algorithm of Spendley et al.,
      \item the variable shape algorithm of Nelder and Mead,
      \item Box's 'complex' algorithm managing bounds and nonlinear inequality
        constraints based on arbitrary number of vertices in the simplex.
      }
    \item Manage various simplex initializations: \itemize{
      \item initial simplex given by user,
      \item initial simplex computed with a length and along the coordinate
        axes,
      \item initial regular simplex computed with formula of Spendley et al.,
      \item initial simplex computed by a small perturbation around the initial
        guess point.
      }
    \item Manage cost function: \itemize{
      \item optional additional argument,
      \item direct communication of the task to perform: cost function or
        inequality constraints.
      }
    \item Manage various termination criteria, including maximum number of
      iterations, tolerance on function value (relative or absolute):\itemize{
      \item tolerance on x (relative or absolute),
      \item tolerance on standard deviation of function value (original
        termination criteria in Box 1965),
      \item maximum number of evaluations of cost function,
      \item absolute or relative simplex size.
      }
    \item Manage the history of the convergence, including: \itemize{
      \item history of function values,
      \item history of optimum point,
      \item history of simplices,
      \item history of termination criteria.
      }
    \item Provide a plot command which allows to graphically see the history
      of the simplices toward the optimum (Not yet implemented).
    \item Provide query features for the status of the optimization process:
      number of iterations, number of function evaluations, status of execution,
      function value at initial point, function value at optimal point, etc...
    \item Kelley restart based on simplex gradient.
    \item O'Neill restart based on factorial search around optimum.
  }
}
\details{
  \tabular{ll}{
  Package: \tab neldermead\cr
  Type: \tab Package\cr
  Version: \tab 1.0-13\cr
  Date: \tab 2026-01-25\cr
  License: \tab CeCILL-2\cr
  LazyLoad: \tab yes
  }
  See \code{vignette('neldermead',package='neldermead')} for more information.
}
\author{
  Author of Scilab neldermead module: Michael Baudin (INRIA - Digiteo)

  Author of R adaptation: Sebastien Bihorel (\email{sb.pmlab@gmail.com})
}
\references{
'Sequential Application of Simplex Designs in Optimisation and Evolutionary
Operation', Spendley, W. and Hext, G. R. and Himsworth, F. R., American
Statistical Association and American Society for Quality, 1962

'A Simplex Method for Function Minimization', Nelder, J. A. and Mead, R., The
Computer Journal, 1965

'A New Method of Constrained Optimization and a Comparison With Other Methods',
M. J. Box, The Computer Journal 1965 8(1):42-52, 1965 by British Computer
Society

'Discussion and correspondence: modification of the complex method of
constrained optimization', J. A. Guin, The Computer Journal, 1968

'Detection and Remediation of Stagnation in the Nelder--Mead Algorithm Using a
Sufficient Decrease Condition', Kelley C. T., SIAM J. on Optimization, 1999

'Iterative Methods for Optimization', C. T. Kelley, SIAM Frontiers in Applied
Mathematics, 1999

'Algorithm AS47 - Function minimization using a simplex procedure', O'Neill, R.,
Applied Statistics, 1971
}
\keyword{ package }
\seealso{
  \code{\link[optimbase]{optimbase}}
  \code{\link[optimsimplex]{optimsimplex}}
}
