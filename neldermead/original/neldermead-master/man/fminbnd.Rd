\name{fminbnd}
\alias{fminbnd}
\title{
  Computation of the constrained minimimum of given function with the
  Nelder-Mead algorithm.
}
\description{
  
  EXPERIMENTAL.
  
  This function searches for the constrained minimum of a given cost function.
  The provided algorithm is a direct search algorithm, i.e. an algorithm which
  does not use the derivative of the cost function. It is based on the update of
  a simplex, which is a set of k>=n+1 vertices, where each vertex is associated
  with one point, which coordinates are constrained within user-defined 
  boundaries, and with one function value. This algorithm corresponds to a 
  version of the Box algorithm, based on bounds and no non-linear constraints. 
  This function is based on a specialized use of the more general
  \code{neldermead} function bundle. Users who want to have a more flexible
  solution based on direct search algorithms should consider using the
  \code{neldermead} functions instead of the \code{fminbnd} function.
}
\usage{
  fminbnd(fun=NULL, x0=NULL, xmin=NULL, xmax=NULL, options=NULL, verbose=FALSE)
}
\arguments{
  \item{fun}{A cost function return a numeric scalar.}
  \item{x0}{A numerical vector of initial guesses (length n).}
  \item{xmin}{A numerical vector of lower bounds for \code{x0} (length n).}
  \item{xmax}{A numerical vector of upper bounds for \code{x0} (length n).}
  \item{options}{A list of optimization options, which drives the behaviour of
    \code{fminbnd}. These options must be set with the \code{optimset}
    function (see \code{?optimset}) which returns a list with the following
    elements: \describe{
      \item{MaxIter}{The maximum number of iterations. The default is 200 * n.}
      \item{MaxFunEvals}{The maximum number of evaluations of the cost function.
        The default is 200 * n.}
      \item{BoxTolFun}{The absolute tolerance on function value. The default value
        is 1.e-4.}
      \item{TolFun}{The absolute tolerance on function value. The default value
        is 1.e-4.}
      \item{TolX}{The absolute tolerance on simplex size. The default value is
        1.e-4.}
      \item{Display}{The verbose level.}
      \item{OutputFcn}{The output function, or a list of output functions
        called at the end of each iteration. The default value is NULL.}
      \item{PlotFcns}{The plot function, or a list of plotput functions called
        at the end of each iteration. The default value is empty.}
    }
  }
  \item{verbose}{The verbose option, controlling the amount of messages.}
}
\details{
  \strong{Termination criteria}

  In this section, we describe the termination criteria used by \code{fminbnd}.
  The criteria is based on the following variables: \describe{
    \item{boxkount}{the current number of time the tolerance on the cost 
      function was met, and}
    \item{shiftfv}{the absolute value of the difference of function value
      between the highest and lowest vertices.}
  }

  If both \code{shiftfv < options$TolFun} and \code{boxkount < options$nbMatch} 
  conditions are true, then the iterations stop.

  \strong{The initial simplex}
  
  The \code{fminbnd} algorithm uses a special initial simplex, which is an
  heuristic depending on the initial guess. The strategy chosen by
  \code{fminbnd} corresponds to the content of \code{simplex0method}
  element of the neldermead object (set to 'randbounds'). It is applied using 
  the content of the \code{boundsmin} and \code{boundsmin} elements to generate
  a simplex with random vertices within the boundaries defined by the user (ie,
  \code{xmin}, and \code{xmax}). This method is an heuristic which is presented 
  in 'A New Method of Constrained Optimization and a Comparison With Other 
  Methods' by M.J. Box. See in the help of \code{optimsimplex} for more details.

  \strong{The number of iterations}
  
  In this section, we present the default values for the number of iterations in
  \code{fminbnd}.

  The \code{options} input argument is an optional list which can
  contain the \code{MaxIter} field, which stores the maximum number of
  iterations. The default value is 200n, where n is the number of variables.
  The factor 200 has not been chosen by chance, but is the result of experiments
  performed against quadratic functions with increasing space dimension.
  This result is presented in 'Effect of dimensionality on the Nelder-Mead
  simplex method' by Lixing Han and Michael Neumann. This paper is based on
  Lixing Han's  PhD, 'Algorithms in Unconstrained Optimization'. The study is
  based on numerical experiments with a quadratic function where the number of
  terms depends on the dimension of the space (i.e. the number of variables).
  Their study showed that the number of iterations required to reach the
  tolerance criteria is roughly 100n. Most iterations are based on inside
  contractions. Since each step of the Nelder-Mead algorithm only require one
  or two function evaluations, the number of required function evaluations in
  this experiment is also roughly 100n.

  \strong{Output and plot functions}
  
  The \code{optimset} function can be used to configure one or more output and
  plot functions.
  The output or plot function is expected to have the following definition:

  \code{myfun <- function(x , optimValues , state)}

  The input arguments \code{x}, \code{optimValues} and \code{state} are
  described in detail in the \code{optimset} help page. The
  \code{optimValues$procedure} field represents the type of step performed at
  the current iteration and can be equal to one of the following strings:
  \itemize{
    \item '' (the empty string),
    \item 'initial simplex',
    \item 'reflect (Box)'.
  }
}
\value{
  Return a object of class neldermead. Use the \code{neldermead.get} to extract
  the following element from the returned object:\describe{
    \item{xopt}{The vector of n numeric values, minimizing the cost function.}
    \item{fopt}{The minimum value of the cost function.}
    \item{exitflag}{The flag associated with exist status of the algorithm.
      The following values are available:\describe{
        \item{-1}{The maximum number of iterations has been reached.}
        \item{0}{The maximum number of function evaluations has been reached.}
        \item{1}{The tolerance on the simplex size and function value delta has
          been reached. This signifies that the algorithm has converged,
          probably to a solution of the problem.}
      }
    }
    \item{output}{A list which stores detailed information about the exit of the
      algorithm. This list contains the following fields:\describe{
        \item{algorithm}{A string containing the definition of the algorithm
          used, i.e. 'Nelder-Mead simplex direct search'.}
        \item{funcCount}{The number of function evaluations.}
        \item{iterations}{The number of iterations.}
        \item{message}{A string containing a termination message.}
      }
    }
  }
}
\author{
Author: Sebastien Bihorel (\email{sb.pmlab@gmail.com})
}
\seealso{
  \code{\link[neldermead]{optimset}}
}
\examples{
#In the following example, we use the fminbnd function to compute the minimum
#of a quadratic function. We first define the function 'quad', and then use
#the fminbnd function to search the minimum, starting with the initial guess
#(1.2, 1.9) and bounds of (1, 1) and (2, 2). In this particular case, 11 
#iterations are performed with 20 function evaluations
  quad <- function(x){
    y <- x[1]^2 + x[2]^2
  }
  sol <- fminbnd(quad,c(1.2,1.9),c(1,1),c(2,2))
  summary(sol)
}
\keyword{ method }
