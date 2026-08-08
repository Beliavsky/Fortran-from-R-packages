\name{neldermead.set}
\alias{neldermead.set}
\title{
  Neldermead Object Configuration
}
\description{
  Configure the current neldermead object with the given value for the given
  key.
}
\usage{
  neldermead.set(this = NULL, key = NULL, value = NULL)
}
\arguments{
  \item{this}{The current neldermead object.}
  \item{key}{The key to configure. See details for the list of possible keys.}
  \item{value}{The value to assign to the key.}
}
\value{
  An updated neldermead object.
}
\details{
  \code{neldermead.set} sets the content of the \code{key} element of
  the neldermead object \code{this} to \code{value}. If \code{key} is a
  sub-element of \code{this$optbase}, \code{value} is assigned by
  \code{optimbase.set}.
  
  The main available keys are the following:\describe{
    \item{'-verbose'}{Set to 1 to enable verbose logging.}
    \item{'-verbosetermination'}{Set to 1 to enable verbose termination
      logging.}
    \item{'-x0'}{The initial guess, as a n x 1 column vector, where n is the
      number of variables.}
    \item{'-maxfunevals'}{The maximum number of function evaluations. If this
      criteria is triggered during optimization, the status of the optimization
      is set to 'maxfuneval'.}
    \item{'-maxiter'}{The maximum number of iterations. If this criteria is
      triggered during optimization, the status of the optimization is set to
      'maxiter'.}option
    \item{'-tolfunabsolute'}{The absolute tolerance for the function value.}
    \item{'-tolfunrelative'}{The relative tolerance for the function value.}
    \item{'-tolfunmethod'}{The method used for the tolerance on function value
      in the termination criteria. The following values are available: TRUE,
      FALSE. If this criteria is triggered, the status of the optimization is
      set to 'tolf'.}
    \item{'-tolxabsolute'}{The absolute tolerance on x.}
    \item{'-tolxrelative'}{The relative tolerance on x.}
    \item{'-tolxmethod'}{The method used for the tolerance on x in the
      termination criteria. The following values are available: TRUE, FALSE. If
      this criteria is triggered during optimization, the status of the
      optimization is set to 'tolx'.}
    \item{'-function'}{The objective function, which computes the value of the
      cost and the non linear constraints, if any. See
      \code{vignette('neldermead',package='neldermead')} for the details of
      the communication between the optimization system and the cost function.}
    \item{'-costfargument'}{An additionnal argument, passed to the cost
      function.}
    \item{'-outputcommand'}{A command which is called back for output. See
      \code{vignette('neldermead',package='neldermead')} for the details of
      the communication between the optimization system and the output command
      function.}
    \item{'-outputcommandarg'}{An additionnal argument, passed to the output
      command.}option
    \item{'-numberofvariables'}{The number of variables to optimize.}
    \item{'-storehistory'}{Set to TRUE to enable the history storing.}
    \item{'-boundsmin'}{The minimum bounds for the parameters.}
    \item{'-boundsmax'}{The maximum bounds for the parameters.}
    \item{'-nbineqconst'}{The number of inequality constraints.}
    \item{'-method'}{The name of the algorithm to use. The following methods are
      available:\describe{
        \item{'fixed'}{the fixed simplex shape algorithm of Spendley et al. This
          algorithm is for unconstrained problems (i.e. bounds and non linear
          constraints are not taken into account)}
        \item{'variable'}{the variable simplex shape algorithm of Nelder and
          Mead. This algorithm is for unconstrained problems (i.e. bounds and
          non linear constraints are not taken into account)}
        \item{'box'}{Box's complex algorithm. This algorithm takes into
          account bounds and nonlinear inequality constraints.}
        \item{'mine'}{the user-defined algorithm, associated with the
          \code{mymethod} element. See
          \code{vignette('neldermead',package='neldermead')} for details.}
      }
    }
    \item{'-simplex0method'}{The method to use to compute the initial simplex.
      The first vertex in the simplex is always the initial guess associated
      with the \code{x0} element. The following methods are available:
      \describe{
        \item{'given'}{The coordinates associated with the \code{coords0}
          element are used to compute the initial simplex, with arbitrary number
          of vertices. This allows the user to setup the initial simplex by a
          specific method which is not provided by the current package (for
          example with a simplex computed from a design of experiments). This
          allows also to configure the initial simplex so that a specific
          behaviour of the algorithm is to be reproduced (for example the Mac
          Kinnon test case). The given matrix is expected to have nbve rows and
          n columns, where n is the dimension of the problem and nbve is the
          number of vertices.}
        \item{'axes'}{The simplex is computed from the coordinate axes and the
          length associated with the \code{simplex0length} element.}
        \item{'spendley'}{The simplex is computed so that it is regular with
          the length associated with the \code{simplex0length} element (i.e. all
          the edges have the same length).}
        \item{'pfeffer'}{The simplex is computed from an heuristic, in the
          neighborhood of the initial guess. This initial simplex depends on
          the -simplex0deltausual and -simplex0deltazero.}
        \item{'randbounds'}{The simplex is computed from the bounds and a
          random number. This option is available only if bounds are available:
          if bounds are not available, an error is generated. This method is
          usually associated with Box's algorithm. The number of vertices in
          the simplex is taken from the \code{boxnbpoints} element.}
      }
    }
    \item{'-coords0'}{The coordinates of the vertices of the initial simplex. If
      the \code{simplex0method} element is set to 'given', these coordinates are
      used to compute the initial simplex. This matrix is expected to have shape
      nbve x n, where nbve is the number of vertices and n is the number of
      variables.}
    \item{'-simplex0length'}{The length to use when the initial simplex is
      computed with the 'axes' or 'spendley' methods. If the initial simplex is
      computed from 'spendley' method, the length is expected to be a
      scalar value. If the initial simplex is computed from 'axes' method,
      it may be either a scalar value or a vector of values, of length n,
      where n is the number of variables.}
    \item{'-simplex0deltausual'}{The relative delta for non-zero parameters in
      'pfeffer' method.}
    \item{'-simplex0deltazero'}{The absolute delta for non-zero parameters in
      'pfeffer' method.}
    \item{'-rho'}{The reflection coefficient. This parameter is used when the
      \code{method} element is set to 'fixed' or 'variable'.}
    \item{'-chi'}{The expansion coefficient. This parameter is used when the
      \code{method} element is set to 'variable'.}
    \item{'-gamma'}{The contraction coefficient. This parameter is used when the
      \code{method} element is set to 'variable'.}
    \item{'-sigma'}{The shrinkage coefficient. This parameter is used when the
      \code{method} element is set to 'fixed' or 'variable'.}
    \item{'-tolsimplexizemethod'}{Set to FALSE to disable the tolerance on the
      simplex size. If this criteria is triggered,
      the status of the optimization is set to 'tolsize'. When this criteria is
      enabled, the values of the \code{tolsimplexizeabsolute} and
      \code{tolsimplexizerelative} elements are used in the termination
      criteria. The method to compute the size is the 'sigmaplus' method.}
    \item{'-tolsimplexizeabsolute'}{The absolute tolerance on the simplex size.}
    \item{'-tolsimplexizerelative'}{The relative tolerance on the simplex size.}
    \item{'-tolssizedeltafvmethod'}{Set to TRUE to enable the termination
      criteria based on the size of the simplex and the difference of function
      value in the simplex. If this criteria is triggered, the status of the
      optimization is set to 'tolsizedeltafv'. This termination criteria uses
      the values of the \code{tolsimplexizeabsolute} and \code{toldeltafv}
      elements.}option
    \item{'-toldeltafv'}{The absolute tolerance on the difference between the
      highest and the lowest function values.}
    \item{'-tolvarianceflag'}{Set to TRUE to enable the termination criteria
      based on the variance of the function value. If this criteria is
      triggered, the status of the optimization is set to 'tolvariance'.
      This criteria is suggested by Nelder and Mead.}
    \item{'-tolabsolutevariance'}{The absolute tolerance on the variance of the
      function values of the simplex.}
    \item{'-tolrelativevariance'}{The relative tolerance on the variance of the
      function values of the simplex.}
    \item{'-kelleystagnationflag'}{Set to TRUE to enable the termination
      criteria using Kelley's stagnation detection, based on sufficient decrease
      condition. If this criteria is triggered, the status of the optimization
      is set to 'kelleystagnation'.}
    \item{'-kelleynormalizationflag'}{Set to FALSE to disable the normalization
      of the alpha coefficient in Kelley's stagnation detection, i.e. use the
      value of the \code{kelleystagnationalpha0} element as is. Default value is
      TRUE, i.e. the simplex gradient of the initial simplex is takeoptionn into
      account in the stagnation detection.}
    \item{'-kelleystagnationalpha0'}{The parameter used in Kelley's stagnation
      detection.}
    \item{'-restartflag'}{Set to TRUE to enable the automatic restart of the
      algorithm.}
    \item{'-restartdetection'}{The method to detect if the automatic restart
      must be performed. The following methods are available:\describe{
      \item{'oneill'}{The factorial local optimality test by O'Neill is used. If
        the test finds a local point which is better than the computed optimum,
        a restart is performed.}
      \item{'kelley'}{The sufficient decrease condition by O'Neill is used. If
        the test finds that the status of the optimization is
        'kelleystagnation', a restart is performed. This status may be generated
        if the -kelleystagnationflag option is set to TRUE.}
      }
    }
    \item{'-restartmax'}{The maximum number of restarts, when automatic restart
      is enabled via the -restartflag option.}
    \item{'-restarteps'}{The absolute epsilon value used to check for optimality
      in the factorial O'Neill restart detection.}
    \item{'-restartstep'}{The absolute step length used to check for optimality
      in the factorial O'Neill restart detection.}
    \item{'-restartsimplexmethod'}{The method to compute the initial simplex
      after a restart. The following methods are available.\describe{
      \item{'given'}{The coordinates associated with the \code{coords0} element
        are used to compute the initial simplex, with arbitrary number of
        vertices. This allow the user to setup the initial simplex by a specific
        method which is not provided by the current package (for example with
        a simplex computed from a design of experiments). This allows also to
        configure the initial simplex so that a specific behaviour of the
        algorithm is to be reproduced (for example the Mc Kinnon test case).
        The given matrix is expected to have nbve rows and n columns, where n
        is the dimension of the problem and nbve is the number of vertices.}
      \item{'axes'}{The simplex is computed from the coordinate axes and the
        length associated with the -simplex0length option.}
      \item{'spendley'}{The simplex is computed so that it is regular with the
        length associated with the -simplex0length option (i.e. all the edges
        have the same length).}
      \item{'pfeffer'}{The simplex is computed from an heuristic, in the
        neighborhood of the initial guess. This initial simplex depends on the
        -simplex0deltausual and -simplex0deltazero.}
      \item{'randbounds'}{The simplex is computed from the bounds and a random
        number. This option is available only if bounds are available: if
        bounds are not available, an error is generated. This method is usually
        associated with Box's algorithm. The number of vertices in the simplex
        is taken from the -boxnbpoints option.}
      \item{'oriented'}{The simplex is computed so that it is oriented, as
        suggested by Kelley.}
      }
    }
    \item{'-scalingsimplex0'}{The algorithm used to scale the initial simplex
      into the nonlinear constraints. The following two algorithms are provided:
      \describe{
      \item{'tox0'}{scales the vertices toward the initial guess.}
      \item{'tocentroid'}{scales the vertices toward the centroid, as
        recommended by Box.}
      }
      If the centroid happens to be unfeasible, because the constraints are not
      convex, the scaling of the initial simplex toward the centroid may fail.
      Since the initial guess is always feasible, scaling toward the initial
      guess cannot fail.}
    \item{'-boxnbpoints'}{The number of points in the initial simplex, when the
      -simplex0method is set to 'randbounds'. The value of this option is also
      use to update the simplex when a restart is performed and the
      -restartsimplexmethod option is set to 'randbounds'. The default value is
      so that the number of points is twice the number of variables of the
      problem.}
    \item{'-boxineqscaling'}{The scaling coefficient used to scale the trial
      point for function improvement or into the constraints of Box's
      algorithm.}
    \item{'-guinalphamin'}{The minimum value of alpha when scaling the vertices
      of the simplex into nonlinear constraints in Box's algorithm.}
    \item{'-boxreflect'}{The reflection factor in Box's algorithm.}
    \item{'-boxtermination'}{Set to TRUE to enable Box's termination criteria.}
    \item{'-boxtolf'}{The absolute tolerance on difference of function values in
      the simplex, suggested by Box. This tolerance is used if the
      -\code{boxtermination} element is set to TRUE.}
    \item{'-boxnbmatch'}{The number of consecutive match of Box's termination
      criteria.}
    \item{'-boxboundsalpha'}{The parameter used to project the vertices into the
      bounds in Box's algorithm.}
    \item{'-mymethod'}{A user-derined simplex algorithm. See
      \code{vignette('neldermead',package='neldermead')} for details.}
    \item{'-myterminate'}{A user-defined terminate function. See
      \code{vignette('neldermead',package='neldermead')} for details.}
    \item{'-myterminateflag'}{Set to TRUE to enable the user-defined terminate
      function.}
  }
}
\author{
Author of Scilab neldermead module: Michael Baudin (INRIA - Digiteo)

Author of R adaptation: Sebastien Bihorel (\email{sb.pmlab@gmail.com})
}
\seealso{
  \code{\link[neldermead]{neldermead}}
}
\keyword{ method }
