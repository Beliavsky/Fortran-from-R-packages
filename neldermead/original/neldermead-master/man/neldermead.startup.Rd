\name{Secondary search functions}
\alias{neldermead.startup}
\alias{neldermead.log}
\alias{neldermead.scaletox0}
\alias{neldermead.scaletocenter}
\alias{neldermead.termstartup}
\alias{neldermead.outputcmd}
\alias{neldermead.autorestart}
\alias{neldermead.istorestart}
\alias{neldermead.isroneill}
\alias{neldermead.isrkelley}
\alias{neldermead.updatesimp}
\alias{scaleinconstraints}
\alias{neldermead.costf}
\title{
  Secondary functions for neldermead.search
}
\description{
  Utility functions for \code{neldermead.serch} and dependent functions.
}
\usage{
  neldermead.startup(this = NULL)
  neldermead.log(this = NULL, msg = NULL)
  neldermead.scaletox0(this = NULL, simplex0 = NULL)
  neldermead.scaletocenter(this = NULL, simplex0 = NULL, x0 = NULL)
  neldermead.termstartup(this = NULL)
  neldermead.outputcmd(this = NULL, state = NULL, simplex = NULL, step = NULL)
  neldermead.autorestart(this = NULL)
  neldermead.istorestart(this = NULL)
  neldermead.isroneill(this = NULL)
  neldermead.isrkelley(this = this)
  neldermead.updatesimp(this = NULL)
  scaleinconstraints(this = NULL, x = NULL, xref = NULL)
  neldermead.costf(x = NULL, this = NULL)
}
\arguments{
  \item{this}{A neldermead object.}
  \item{msg}{A character string.}
  \item{simplex0}{The initial simplex object.}
  \item{x0}{A column matrix of initial parameters.}
  \item{state}{The state of the algorithm, either 'init', 'done' or 'iter'.}
  \item{simplex}{The current simplex object.}
  \item{step}{The type of step performed during the iteration: 'init', 'done',
    'reflection', 'expansion', 'insidecontraction', 'outsidecontraction',
    'reflectionnext' or 'shrink'.}
  \item{x}{The point estimate to scale.}
  \item{xref}{The reference point estimate.}
}
\details{
\describe{
  \item{\code{neldermead.startup}}{Startup the algorithm. Compute the
    initial simplex, depending on the content of the \code{simplex0method}
    element of the neldermead object ('given', 'axes', 'spendley', 'pfeffer' or
    'randbounds').}
  \item{\code{neldermead.log}}{Print a message to the log file using
    \code{optimbase.log}.}
  \item{\code{neldermead.scaletox0}}{Scale the simplex into the nonlinear
    inequality constraints, if any. Scale toward x0, which is feasible.}
  \item{\code{neldermead.scaletocenter}}{Scale the simplex into the nonlinear
    inequality constraints, if any. Scale to the centroid of the points which
    satisfy the constraints. This is Box's method for scaling. It is unsure,
    since the centroid of the points which satisfy the constraints may not
    be feasible.}
  \item{\code{neldermead.termstartup}}{Initialize Kelley's stagnation detection
    system when normalization is required, by computing kelleyalpha. If the
    simplex gradient is zero, then use alpha0 as alpha.}
  \item{\code{neldermead.outputcmd}}{Call the array of user-defined output
    functions}
  \item{\code{neldermead.autorestart}}{Perform an optimization with automatic
    restart. The loop processes for i = 1 to \code{restartmax} + 1. This is
    because a RE-start is performed after one simulation has been performed,
    hence the 'RE'.}
  \item{\code{neldermead.istorestart}}{Determine if the optimization is to
    restart using \code{neldermead.isroneill} or \code{neldermead.isrkelley}
    depending on the content of the \code{restartdetection} element.}
  \item{\code{neldermead.isroneill}}{Determine if the optimization is to
    restart. Use O'Neill method as a criteria for restart. It is an axis-by-axis
    search for optimality.}
  \item{\code{neldermead.isrkelley}}{Determine if the optimization is to
    restart. Use \code{kelleystagnation} as a criteria for restart.}
  \item{\code{neldermead.updatesimp}}{Update the initial simplex simplex0 for a
    restart.}
  \item{\code{scaleinconstraints}}{Given a point reference to scale and a
    reference point which satisfies the constraints, scale the point
    towards the reference point estimate until it satisfies all the
    constraints.}
  \item{\code{neldermead.costf}}{Call the cost function and return the value.
    This function is given to the simplex function class as a callback.
    Input/Output arguments are swapped w.r.t. \code{optimbase.function}, so that
    it matches the requirements of simplex methods.}
}
}
\value{
\describe{
  \item{\code{neldermead.startup}}{Return an updated neldermead object
    \code{this}.}
  \item{\code{neldermead.log}}{Return the neldermead object \code{this}.}
  \item{\code{neldermead.scaletox0}}{Return an updated simplex.}
  \item{\code{neldermead.scaletocenter}}{Return an updated simplex.}
  \item{\code{neldermead.termstartup}}{Return an updated neldermead object
    \code{this}.}
  \item{\code{neldermead.outputcmd}}{Do not return any data, but execute the
    output function(s).}
  \item{\code{neldermead.autorestart}}{Return an updated neldermead object
    \code{this}.}
  \item{\code{neldermead.istorestart}}{Return a list with the following
    elements:\describe{
      \item{this}{The input neldermead object.}
      \item{istorestart}{Set to TRUE if the optimization is to restart, to FALSE
        otherwise.}
    }
  }
  \item{\code{neldermead.isroneill}}{Return a list with the following elements:
    \describe{
      \item{this}{The input neldermead object.}
      \item{istorestart}{Set to TRUE if the optimization is to restart, to FALSE
        otherwise.}
    }
  }
  \item{\code{neldermead.isrkelley}}{Return a list with the following elements:
    \describe{
      \item{this}{The input neldermead object.}
      \item{istorestart}{Set to TRUE if the optimization is to restart, to FALSE
        otherwise.}
    }
  }
  \item{\code{neldermead.updatesimp}}{Return an updated neldermead object
    \code{this}.}
  \item{\code{scaleinconstraints}}{Return a list with the following elements:
    \describe{
      \item{this}{The updated neldermead object.}
      \item{isscaled}{TRUE if the procedure has succeeded before
        \code{boxnbnlloops}, FALSE if it has failed.}
      \item{p}{The scaled parameters.}
    }
  }
  \item{\code{neldermead.costf}}{Return a list with the following elements:
    \describe{
      \item{f}{The value of the cost function.}
      \item{this}{The updated neldermead object.}
    }
  }
}
}

\author{
Author of Scilab neldermead module: Michael Baudin (INRIA - Digiteo)

Author of R adaptation: Sebastien Bihorel (\email{sb.pmlab@gmail.com})
}
\keyword{ method }
