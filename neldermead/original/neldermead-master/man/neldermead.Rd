\name{neldermead}
\alias{neldermead}
\alias{print.neldermead}
\alias{summary.neldermead}
\alias{is.neldermead}
\title{
  Neldermead objects
}
\description{
  These functions support the S3 class 'neldermead' and are intended to 
  either create objects of this class or check if an object is of this class.
}
\usage{
  neldermead(optbase, method, simplex0, simplex0method, 
    simplex0length, simplexsize0, simplexopt, historysimplex, coords0, rho, chi, 
    gamma, sigma, tolfstdeviation, tolfstdeviationmethod, tolsimplexizeabsolute, 
    tolsimplexizerelative, tolsimplexizemethod, toldeltafv, tolssizedeltafvmethod,
    simplex0deltausual, simplex0deltazero, restartsimplexmethod, restartmax, 
    restarteps, restartstep, restartnb, restartflag, restartdetection, 
    kelleystagnationflag, kelleynormalizationflag, kelleystagnationalpha0, 
    kelleyalpha, startupflag, boxnbpoints, boxnbpointseff, boxineqscaling, 
    checkcostfunction, scalingsimplex0, guinalphamin, boxboundsalpha, 
    boxtermination, boxtolf, boxnbmatch, boxkount, boxreflect, tolvarianceflag, 
    tolabsolutevariance, tolrelativevariance, variancesimplex0, mymethod, 
    myterminate, myterminateflag, greedy, output, exitflag)
    
  \method{print}{neldermead}(x,verbose,...)
  
  \method{summary}{neldermead}(object,showhistory,...)
  
  \method{is}{neldermead}(x=NULL)
  
}
\arguments{
  \item{optbase}{An object of class 'optimbase', i.e. a list created by 
    \code{optimbase()} and containing the following elements: \describe{
      \item{verbose}{The verbose option, controlling the amount of messages.}
      \item{x0}{The initial guess.}
      \item{fx0}{The value of the function for the initial guess.}
      \item{xopt}{The optimum parameter.}
      \item{fopt}{The optimum function value.}
      \item{tolfunabsolute}{The absolute tolerance on function value.}
      \item{tolfunrelative}{The relative tolerance on function value.}
      \item{tolfunmethod}{Logical flag for the tolerance on function value in
        the termination criteria. This criteria is suitable for functions which
        minimum is associated with a function value equal to 0.}
      \item{tolxabsolute}{The absolute tolerance on x.}
      \item{tolxrelative}{The relative tolerance on x.}
      \item{tolxmethod}{Possible values: FALSE, TRUE.}
      \item{funevals}{The number of function evaluations.}
      \item{maxfunevals}{The maximum number of function evaluations.}
      \item{maxiter}{The maximum number of iterations.}
      \item{iterations}{The number of iterations.}
      \item{fun}{The cost function.}
      \item{status}{The status of the optimization.}
      \item{historyfopt}{The vector to store the history for fopt. The values of
        the cost function will be stored at each iteration in a new element, so
        the length of \code{historyfopt} at the end of the optimization should
        be the number of iterations.}
      \item{historyxopt}{The list to store the history for xopt. The vectors of 
        estimates will be stored on separated levels of the list, so the length
        of \code{historyfopt} at the end of the optimization should be the
        number of iterations.}
      \item{verbosetermination}{The verbose option for termination criteria.}
      \item{outputcommand}{The command called back for output.}
      \item{outputcommandarg}{The outputcommand argument is initialized as a
        string. If the user configure this element, it is expected that a matrix
        of values or a list is passed so that the argument is appended to the
        name of the function.}
      \item{numberofvariables}{The number of variables to optimize.}
      \item{storehistory}{The flag which enables/disables the storing of the
        history.}
      \item{costfargument}{The costf argument is initialized as a string. If
        the user configure this element, it is expected that a matrix of values
        or a list is passed so that the argument is appended to the name of the
        function.}
      \item{boundsmin}{Minimum bounds for the parameters.}
      \item{boundsmax}{Maximum bounds for the parameters.}
      \item{nbineqconst}{The number of nonlinear inequality constraints.}
      \item{logfile}{The name of the log file.}
      \item{logfilehandle}{The handle for the log file.}
      \item{logstartup}{Set to TRUE when the logging is started up.}
      \item{withderivatives}{Set to TRUE when the method uses derivatives.}
    }
  }
  \item{method}{The name of the algorithm to use.}
  \item{simplex0}{An object of class 'simplex', i.e. a list created by
    \code{optimsimplex()}, and containing the following elements:\describe{
      \item{verbose}{The verbose option, controlling the amount of messages.}
      \item{x}{The coordinates of the vertices, with size nbve x n.}
      \item{n}{The dimension of the space.}
      \item{fv}{The function values, with size nbve x 1.}
      \item{nbve}{The number of vertices.}
    }
  }
  \item{simplex0method}{The method to use to compute the initial simplex.}
  \item{simplex0length}{The length to use when the initial simplex is computed
    with the 'axes' or 'spendley' methods.}
  \item{rho}{The reflection coefficient. This parameter is used when the
    \code{method} element is set to 'fixed' or 'variable'.}
  \item{chi}{The expansion coefficient. This parameter is used when the
    \code{method} element is set to 'variable'.}
  \item{gamma}{The contraction coefficient. This parameter is used when the
    \code{method} element is set to 'variable'.}
  \item{sigma}{The shrinkage coefficient. This parameter is used when the
    \code{method} element is set to 'fixed' or 'variable'.}
  \item{tolfstdeviation}{The tolerance for the standard deviation.}
  \item{tolfstdeviationmethod}{Set to FALSE.}
  \item{tolsimplexizeabsolute}{The absolute tolerance on the simplex size.}
  \item{tolsimplexizerelative}{The relative tolerance on the simplex size.}
  \item{tolsimplexizemethod}{Logical flag to enable/disable the tolerance on
    the simplex size. When this criteria is enabled, the values of the
    \code{tolsimplexizeabsolute} and \code{tolsimplexizerelative} elements are
    used in the termination criteria. The method to compute the size is the
    'sigmaplus' method.}
  \item{simplexsize0}{Initial size of the simplex, for the tolerance on the
    simplex size.}
  \item{toldeltafv}{The absolute tolerance on the difference between the
    highest and the lowest function values.}
  \item{tolssizedeltafvmethod}{Logical flag to enable/disable the termination
    criteria based on the size of the simplex and the difference of function
    value in the simplex. If this criteria is triggered, the status of the
    optimization is set to 'tolsizedeltafv'. This termination criteria uses
    the values of the \code{tolsimplexizeabsolute} and \code{toldeltafv}
    elements. This criteria is identical to Scilab's \code{fminsearch}.}
  \item{historysimplex}{The list to store the history for simplex. The simplex
    will be stored on a new level of the list at each iteration, so the length 
    of \code{historyfopt} at the end of the optimization should be the number
    of iterations.}
  \item{coords0}{The coordinates of the vertices of the initial simplex. If
    the \code{simplex0method} element is set to 'given', these coordinates are
    used to compute the initial simplex. This matrix is expected to have
    shape nbve x n where nbve is the number of vertices and n is the number of
    variables.}
  \item{simplex0deltausual}{The relative delta for non-zero parameters in
    'pfeffer' method.}
  \item{simplex0deltazero}{The absolute delta for non-zero parameters in
    'pfeffer' method.}
  \item{simplexopt}{The optimum simplex, after one optimization process.}
  \item{restartsimplexmethod}{The method to compute the initial simplex after
    a restart.}
  \item{restartmax}{The maximum number of restarts, when automatic restart is
    enabled via the \code{restartflag} element.}
  \item{restarteps}{The absolute epsilon value used to check for optimality in
    the factorial O'Neill restart detection.}
  \item{restartstep}{The absolute step length used to check for optimality in
    the factorial O'Neill restart detection.}
  \item{kelleystagnationflag}{Logical flag to enable/disable the termination
    criteria using Kelley's stagnation detection, based on sufficient decrease
    condition. If this criteria is triggered, the status of the optimization
    is set to 'kelleystagnation'.},
  \item{kelleynormalizationflag}{Logical flag to enable/disable the
    normalization of the alpha coefficient in Kelley's stagnation detection,
    i.e. use the value of the \code{kelleystagnationalpha0} element as is.}
  \item{kelleystagnationalpha0}{The parameter used in Kelley's stagnation
    detection.}
  \item{kelleyalpha}{The current value of Kelley's alpha, after normalization,
    if required.}
  \item{restartnb}{Number of restarts performed.}
  \item{restartflag}{Logical flag to enable/disable the automatic restart of
    the algorithm.}
  \item{restartdetection}{The method to detect if the automatic restart must
    be performed.}
  \item{startupflag}{Set to TRUE when the startup has been performed.}
  \item{boxnbpoints}{The number of points in the initial simplex, when the
    \code{simplex0method} is set to 'randbounds'. The value of this element is
    also use to update the simplex when a restart is performed and the
    \code{restartsimplexmethod} element is set to 'randbounds'. The default
    value is so that the number of points is twice the number of variables of
    the problem.}
  \item{boxnbpointseff}{The effective number of points required in the simplex
    for Box's algorithm.}
  \item{boxineqscaling}{The scaling coefficient used to scale the trial point
    for function improvement or into the constraints of Box's algorithm.}
  \item{checkcostfunction}{Logical flag to enable/disable the checking of the
    connection of the cost function.}
  \item{scalingsimplex0}{The algorithm used to scale the initial simplex into
    the nonlinear constraints. The following two algorithms are provided:
    \describe{
    \item{'tox0'}{scales the vertices toward the initial guess.}
    \item{'tocentroid'}{scales the vertices toward the centroid, as
      recommended by Box.}
    }
    If the centroid happens to be unfeasible, because the constraints are not
    convex, the scaling of the initial simplex toward the centroid may fail.
    Since the initial guess is always feasible, scaling toward the initial
    guess cannot fail.}
  \item{guinalphamin}{The minimum value of alpha when scaling the vertices of
    the simplex into nonlinear constraints in Box's algorithm.}
  \item{boxboundsalpha}{The parameter used to project the vertices into the
    bounds in Box's algorithm.}
  \item{boxtermination}{Logical flag to enable/disable Box's termination
    criteria.}
  \item{boxtolf}{The absolute tolerance on difference of function values in
    the simplex, suggested by Box. This tolerance is used if the
    \code{boxtermination} element is set to TRUE.}
  \item{boxnbmatch}{The number of consecutive match of Box's termination
    criteria.}
  \item{boxkount}{Current number of consecutive match.}
  \item{boxreflect}{The reflection factor in Box's algorithm.}
  \item{tolvarianceflag}{Logical flag to enable/disable the termination
    criteria based on the variance of the function value. If this criteria is
    triggered, the status of the optimization is set to 'tolvariance'. This
    criteria is suggested by Nelder and Mead.}
  \item{tolabsolutevariance}{The absolute tolerance on the variance of the
    function values of the simplex.}
  \item{tolrelativevariance}{The relative tolerance on the variance of the
    function values of the simplex.}
  \item{variancesimplex0}{Relative tolerance on variance.}
  \item{mymethod}{A user-derined simplex algorithm.}
  \item{myterminate}{A user-defined terminate function.}
  \item{myterminateflag}{Logical flag to enable/disable the user-defined
    terminate function.}
  \item{greedy}{Logical flag to enable/disable greedy Nelder-Mead.}
  \item{output}{The command to call back for user-defined output of specialized 
    function.}
  \item{exitflag}{Logical flag to enable/disable the user-defined output of
    specialized function.}
  \item{x}{An object of class 'neldermead'.}
  \item{verbose}{A logical flag, controlling the amount of data printed.}
  \item{...}{optional arguments to 'print' or 'plot' methods.}
  \item{object}{An object of class 'neldermead'.}
  \item{showhistory}{Optional logical flag, to define whether optimization 
    history must be summarized or not.}
}
\value{
  The \code{neldermead} function returns a new object of class 'neldermead', 
  with the following default content:\describe{
    \item{optbase}{An object of class 'optimbase' with the following default 
      content: \describe{
        \item{verbose}{Default is FALSE.}
        \item{x0}{Default is NULL.}
        \item{fx0}{Default is NULL.}
        \item{xopt}{Default is 0.}
        \item{fopt}{Default is 0.}
        \item{tolfunabsolute}{Default is 0.}
        \item{tolfunrelative}{Default is .Machine$double.eps.}
        \item{tolfunmethod}{Default is FALSE.}
        \item{tolxabsolute}{Default is 0.}
        \item{tolxrelative}{Default is .Machine$double.eps.}
        \item{tolxmethod}{Default is TRUE.}
        \item{funevals}{Default is 0.}
        \item{maxfunevals}{Default is 100.}
        \item{maxiter}{Default is 100.}
        \item{iterations}{Default is 0.}
        \item{fun}{Default is ''.}
        \item{status}{Default is ''.}
        \item{historyfopt}{Default is NULL.}
        \item{historyxopt}{Default is NULL.}
        \item{verbosetermination}{Default is FALSE.}
        \item{outputcommand}{Default is ''.}
        \item{outputcommandarg}{Default is ''. If the user configures this 
          element, it is expected to be an object of class 'optimbase.outputargs' 
          or will be coerced to an object of class 'optimbase.outputargs'.}
        \item{numberofvariables}{Default is 0.}
        \item{storehistory}{Default is FALSE.}
        \item{costfargument}{Default is ''. If the user configures this element, 
          it is expected to be an object of class 'optimbase.functionargs' or 
          will be coerced to an object of class 'optimbase.functionargs'.}
        \item{boundsmin}{Default is NULL.}
        \item{boundsmax}{Default is NULL.}
        \item{nbineqconst}{Default is 0.}
        \item{logfile}{Default is ''.}
        \item{logfilehandle}{Default is 0.}
        \item{logstartup}{Default is FALSE.}
        \item{withderivatives}{Default is FALSE.}
      }
    }
    \item{method}{Default is 'variable'.}
    \item{simplex0}{Default is an object of class 'simplex', with the following 
      content:\describe{
        \item{verbose}{Default is 0.}
        \item{x}{Default is NULL.}
        \item{n}{Default is 0.}
        \item{fv}{Default is NULL.}
        \item{nbve}{Default is 0.}
      }
    }
    \item{simplex0method}{Default is 'axes'.}
    \item{simplex0length}{Default is 1.}
    \item{rho}{Default is 1.}
    \item{chi}{Default is 2.}
    \item{gamma}{Default is 0.5.}
    \item{sigma}{Default is 0.5.}
    \item{tolfstdeviation}{Default is 0.}
    \item{tolfstdeviationmethod}{Default is FALSE.}
    \item{tolsimplexizeabsolute}{Default is 0.}
    \item{tolsimplexizerelative}{Default is .Machine$double.eps.}
    \item{tolsimplexizemethod}{Default is FALSE.}
    \item{simplexsize0}{Default is 0.}
    \item{toldeltafv}{Default is .Machine$double.eps.}
    \item{tolssizedeltafvmethod}{Default is FALSE.}
    \item{historysimplex}{Default is NULL.}
    \item{coords0}{Default is NULL.}
    \item{simplex0deltausual}{Default is 0.05.}
    \item{simplex0deltazero}{Default is 0.0075.}
    \item{simplexopt}{Default is NULL.}
    \item{restartsimplexmethod}{Default is 'oriented'.}
    \item{restartmax}{Default is 3.}
    \item{restarteps}{Default is .Machine$double.eps.}
    \item{restartstep}{Default is 1.}
    \item{kelleystagnationflag}{Default is FALSE.},
    \item{kelleynormalizationflag}{Default is TRUE, i.e. the simplex gradient of 
      the initial simplex is taken into account in the stagnation detection.}
    \item{kelleystagnationalpha0}{Default is 1.e-4.}
    \item{kelleyalpha}{Default is 1.e-4.}
    \item{restartnb}{Default is 0.}
    \item{restartflag}{Default is FALSE.}
    \item{restartdetection}{Default is 'oneill'.}
    \item{startupflag}{Default is FALSE.}
    \item{boxnbpoints}{Default is '2n'.}
    \item{boxnbpointseff}{Default is 0.}
    \item{boxineqscaling}{Default is 0.}
    \item{checkcostfunction}{Default is TRUE.}
    \item{scalingsimplex0}{Default is 'tox0'.}
    \item{guinalphamin}{Default is 1.e-6.}
    \item{boxtermination}{Default is FALSE.}
    \item{boxtolf}{Default is 1.e-5.}
    \item{boxnbmatch}{Default is 5.}
    \item{boxkount}{Default is 0.}
    \item{boxreflect}{Default is 1.3.}
    \item{tolvarianceflag}{Default is FALSE.}
    \item{tolabsolutevariance}{Default is 0.}
    \item{tolrelativevariance}{Default is .Machine$double.eps.}
    \item{variancesimplex0}{Default is .Machine$double.eps.}
    \item{mymethod}{Default is NULL.}
    \item{myterminate}{Default is NULL.}
    \item{myterminateflag}{Default is FALSE.}
    \item{greedy}{Default is FALSE.}
    \item{output}{Default is list().}
    \item{exitflag}{Default is FALSE.}
  }
}
\author{
Author of Scilab neldermead module: Michael Baudin (INRIA - Digiteo)

Author of R adaptation: Sebastien Bihorel (\email{sb.pmlab@gmail.com})
}
\seealso{
  \code{\link[optimbase]{optimbase}},
  \code{\link[optimsimplex]{optimsimplex}}
}
\keyword{ method }
