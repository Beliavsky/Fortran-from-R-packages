\name{neldermead.algo}
\alias{neldermead.algo}
\alias{neldermead.fixed}
\alias{neldermead.variable}
\alias{neldermead.box}
\alias{neldermead.storehistory}
\alias{neldermead.termination}
\alias{neldermead.interpolate}
\alias{boxlinesearch}
\title{
  Nelder-Mead Algorithm
}
\description{
  \code{neldermead.algo} performs an optimization without restart using the
  method associated with the \code{method} element of the neldermead object;
  \code{neldermead.fixed}, \code{neldermead.variable}, \code{neldermead.box},
  \code{boxlinesearch}, \code{neldermead.storehistory},
  \code{neldermead.termination}, and \code{neldermead.interpolate} are utility
  functions for \code{neldermead.algo}.
}
\usage{
  neldermead.algo(this = NULL)
  neldermead.fixed(this = NULL)
  neldermead.variable(this = NULL)
  neldermead.box(this = this)
  boxlinesearch(this = NULL, n = NULL, xbar = NULL, xhigh = NULL, fhigh = NULL,
                rho = NULL)
  neldermead.storehistory(this = NULL, n = NULL, fopt = NULL, xopt = NULL,
                          fv = NULL, xcoords = NULL)
  neldermead.termination(this = NULL, fvinitial = NULL, oldfvmean = NULL,
                         newfvmean = NULL, previousxopt = NULL,
                         currentxopt = NULL, simplex = NULL)
  neldermead.interpolate(x1 = NULL, x2 = NULL, fac = NULL)
}
\arguments{
  \item{this}{A neldermead object.}
  \item{n}{Number of variables.}
  \item{xbar}{The centroid.}
  \item{xhigh}{The high point.}
  \item{fhigh}{The value of the cost function at \code{xhigh}.}
  \item{rho}{The reflection factor.}
  \item{fopt}{The current value of the function at the current optimum point
    estimate.}
  \item{xopt}{The current optimum point estimate.}
  \item{fv}{The function values, with size nbve x 1.}
  \item{xcoords}{Matrix of size n x n+1, coordinates of the n+1 vertices}
  \item{fvinitial}{The initial cost function value.}
  \item{oldfvmean}{The old cost function value average on the simplex.}
  \item{newfvmean}{The new cost function value average on the simplex.}
  \item{previousxopt}{The previous point estimate.}
  \item{currentxopt}{The current point estimate.}
  \item{simplex}{The simplex. The best point estimate in the simplex is expected
    to be stored at 1, while the worst point estimate in the simplex is expected
    to be stored at n+1.}
  \item{x1}{The first reference point estimate to perform the
    interpolation.}
  \item{x2}{The second reference point estimate to perform the
    interpolation.}  
  \item{fac}{A factor to perform the interpolation.}
}
\details{
\describe{
  \item{\code{neldermead.fixed}}{The simplex algorithm with fixed size simplex.
    We implement the following 'rules' of the method of Spendley et al.
    \itemize{
    \item Rule 1 is strictly applied, but the reflection is done by reflection
      of the high point, since we minimize a function instead of maximizing it,
      like Spendley.
    \item Rule 2 is NOT implemented, as we expect that the function evaluation
      is not subject to errors.
    \item Rule 3 is applied, i.e. reflection with respect to next to high point.
       A shrink step is included, with shrinkage factor sigma.
    }

    Rule 1. Ascertain the lowest reading y, of yi ... Yk+1
    Complete a new simplex Sp by excluding the point Vp corresponding to
    y, and replacing it by V* defined as above.

    Rule 2. If a result has occurred in (k + 1) successive simplexes, and is not
    then eliminated by application of Rule 1, do not move in the direction
    indicated by Rule 1, or at all, but discard the result and replace it by a
    new observation at the same point.
    
    Rule 3. If y is the lowest reading in So , and if the next observation made,
    y* , is the lowest reading in the new simplex S , do not apply Rule 1 and
    return to So from Sp . Move out of S, by rejecting the second lowest reading
    (which is also the second lowest reading in So).
  }
  \item{\code{neldermead.variable}}{The original Nelder-Mead algorithm, with
    variable-size simplex.}
  \item{\code{neldermead.box}}{The Nelder-Mead algorithm, with variable-size
    simplex and modifications by Box for bounds and inequality constraints.}
  \item{\code{boxlinesearch}}{Called by \code{neldermead.box}, i.e. Box's
    method. Perform a line search from xbar, on the line (xhigh,xbar). The
    reflected point estimate satisfies the following constraints: \itemize{
      \item fr < fhigh
      \item xr satisfies the bounds constraints
      \item xr satisfies the nonlinear positive inequality constraints
      \item xr satisfies the linear positive inequality constraints
    }
    The method is based on projection and scaling toward the centroid.}
  \item{\code{neldermead.storehistory}}{Store the optimization history into the
    neldermead object.}
  \item{\code{neldermead.termination}}{Determine if the algorithm must continue
    or terminate. The function uses the cost function average in the simplex
    instead of the best cost function value. This is because the function
    average changes at each iteration. Instead, the best function value has a
    step-by-step evolution and may not change between two successive iterations,
    leading to a stop of the algorithm.}
  \item{\code{neldermead.interpolate}}{Compute the point estimate xi as an
    interpolation between x1 and x2, as follows:  xi = (1+fac)x1 - fac*x2}
}
}
\value{
\describe{
  \item{\code{neldermead.fixed}, \code{neldermead.variable}, and
    \code{neldermead.box}}{Return the updated neldermead object, containing the
    optimum point estimate.}
  \item{\code{boxlinesearch}}{Return a list with the following elements:
    \describe{
    \item{this}{The updated neldermead object.}
    \item{status}{TRUE if the search is successful, FALSE otherwise.}
    \item{xr}{The reflected point estimate.}
    \item{fr}{The value of the cost function at \code{xr}.}
    }
  }
  \item{\code{neldermead.storehistory}}{Return the updated neldermead object.}
  \item{\code{neldermead.termination}}{Return a list with the following
    elements: \describe{
      \item{this}{The updated neldermead object}
      \item{terminate}{TRUE if the algorithm terminates, FALSE if the algorithm
        must continue.}
      \item{status}{The termination status: 'continue', 'maxiter', 'maxfuneval',
        'tolf', 'tolx', 'tolsize', 'tolsizedeltafv', 'kelleystagnation',
        'tolboxf', 'tolvariance' or the user-defined termination status.}
    }
  }
  \item{\code{neldermead.interpolate}}{Return a new point estimate, i.e. a
    column vector.}
}
}
\author{
  Author of Scilab neldermead module: Michael Baudin (INRIA - Digiteo)

  Author of R adaptation: Sebastien Bihorel (\email{sb.pmlab@gmail.com})
}

\keyword{ method }
