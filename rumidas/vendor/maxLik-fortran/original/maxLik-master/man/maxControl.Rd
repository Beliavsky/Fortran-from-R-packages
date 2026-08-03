\name{MaxControl-class}
\Rdversion{1.1}
\docType{class}
\alias{MaxControl-class}
\alias{maxControl}
\alias{maxControl,MaxControl-method}
\alias{maxControl,missing-method}
\alias{maxControl,maxim-method}
\alias{show,MaxControl-method}

\title{Class \code{"MaxControl"}}
\description{
  This is the structure that holds the optimization control options.
  The corresponding constructors take
  the parameters, perform consistency checks, and return the
  control structure.  Alternatively, it overwrites the supplied
  parameters in an existing \code{MaxControl} structure.  There is also
  a method to extract the control structure from the estimated
  \sQuote{maxim}-objects.  
}
\section{Slots}{
  The default values and definition of the slots:
  \describe{
    \item{tol}{1e-8, stopping condition 
      for \code{\link{maxNR}} and related optimizers.
      Stop if the absolute difference
      between successive iterations is less than \code{tol}, returns
      code 2.}
    \item{reltol}{sqrt(.Machine$double.eps), relative convergence
      tolerance (used by \code{\link{maxNR}} related optimizers, and
      \code{\link{optim}}-based optimizers.
      The algorithm stops if
      it iteration increases the value by less than a factor of
      \code{reltol*(abs(val) + reltol)}.
      Returns code 2.}
    \item{gradtol}{1e-6, stopping condition
      for \code{\link{maxNR}} and related optimizers.
      Stops if norm of the gradient is
      less than \code{gradtol}, returns code 1.}
    \item{steptol}{1e-10, stopping/error condition
      for \code{\link{maxNR}} and related optimizers.
      If \code{qac == "stephalving"} and the quadratic
      approximation leads to a worse, instead of a better value, or to
      \code{NA}, the step length
      is halved and a new attempt is made.  If necessary, this procedure is repeated
      until \code{step < steptol}, thereafter code 3 is returned.}
    %
    \item{lambdatol}{1e-6, (for \code{\link{maxNR}} related
      optimizers)
      controls whether Hessian is treated as negative
      definite.  If the
      largest of the eigenvalues of the Hessian is larger than
      \code{-lambdatol} (Hessian is not negative definite),
      a suitable diagonal matrix is subtracted from the
      Hessian (quadratic hill-climbing) in order to enforce negative
      definiteness.}
    %
    \item{qac}{"stephalving", character, Qadratic Approximation
      Correction for \code{\link{maxNR}} related optimizers.  When the new
      guess is worse than the initial one, program attempts to correct it:
      \code{"stephalving"} decreases the
      step but keeps the direction.
      \code{"marquardt"} uses
      \cite{Marquardt (1963)} method by decreasing the step length while also
      moving closer to the pure gradient direction.  It may be faster and
      more robust choice in areas where quadratic approximation behaves poorly.}
    \item{qrtol}{1e-10, QR-decomposition tolerance
      for Hessian inversion in \code{\link{maxNR}} related optimizers.
    }
    \item{marquardt_lambda0}{0.01, a positive numeric, initial correction term
      for \cite{Marquardt (1963)} correction in
      \code{\link{maxNR}}-related optimizers}
    \item{marquardt_lambdaStep}{2, how much the \cite{Marquardt (1963)}
      correction is decreased/increased at
      successful/unsuccesful step
      for \code{\link{maxNR}} related optimizers}
    \item{marquardt_maxLambda}{1e12, maximum allowed correction term
      for \code{\link{maxNR}} related optimizers.
      If exceeded, the
      algorithm exits with return code 3.}
    %
    \item{nm_alpha}{1, Nelder-Mead simplex method reflection
      factor (see Nelder & Mead, 1965)}
    \item{nm_beta}{0.5, Nelder-Mead contraction factor}
    \item{nm_gamma}{2, Nelder-Mead expansion factor}
    % SANN
    \item{sann_cand}{\code{NULL} or a function for \code{"SANN"} algorithm
      to generate a new candidate point;
      if \code{NULL}, Gaussian Markov kernel is used
      (see argument \code{gr} of \code{\link{optim}}).}
    \item{sann_temp}{10, starting temperature
      for the \dQuote{SANN} cooling schedule.  See \code{\link{optim}}.}
    \item{sann_tmax}{10, number of function evaluations at each temperature for
      the \dQuote{SANN} optimizer.  See \code{\link{optim}}.}
    \item{sann_randomSeed}{123, integer to seed random numbers to
      ensure replicability of \dQuote{SANN} optimization and preserve
      \code{R} random numbers.  Use
      options like \code{SANN_randomSeed=Sys.time()} or
      \code{SANN_randomeSeed=sample(1000,1)} if you want stochastic results.
    }
    % SG general
    General options for stochastic gradient methods:
    \item{SG_learningRate}{0.1, learning rate, numeric}
    \item{SG_batchSize}{\code{NULL}, batch size for Stochastic Gradient Ascent.  A
      positive integer, or \code{NULL} for full-batch gradent ascent.}
    \item{SG_clip}{\code{NULL}, gradient clipping threshold.  This is
      the max allowed squared Euclidean norm of the gradient.  If the
      actual norm of the gradient exceeds (square root of) this
      threshold, the gradient will be scaled back accordingly while
      preserving its direction.  \code{NULL} means no clipping.
    }
    \item{SG_patience}{\code{NULL}, or integer.  Stopping condition: if
      the objective function is worse than its largest value so far this
      many times, the algorithm stops, and returns not the last
      parameter value but the one that
      gave the best results so far.  This is mostly useful if gradient
      is computed on training data and the
      objective function on validation data.
    }
    \item{SG_patienceStep}{1L, integer.  After how many epochs to check
      the patience value.  1 means to check (and hence to compute the
      objective function) at each epoch.
    }
    % Stochastic Gradient Ascent
    Options for SGA:
    \item{SGA_momentum}{0, numeric momentum parameter for SGA.  Must lie
      in interval \eqn{[0,1]}{[0,1]}.
    }
    % Adam
    Options for Adam:
    \item{Adam_momentum1}{0.9, numeric in \eqn{[0,1]}{[0,1]}, the first moment momentum}
    \item{Adam_momentum2}{0.999, numeric in \eqn{[0,1]}{[0,1]}, the second moment momentum}
    % general
    General options:
    \item{iterlim}{150, stopping condition (the default differs for
      different methods).  Stop if more than \code{iterlim}
      iterations performed.  Note that \sQuote{iteration} may mean
      different things for different optimizers.}
    \item{max.rows}{20, maximum number of matrix rows to be printed when
      requesting verbosity in the optimizers.
    }
    \item{max.cols}{7, maximum number of columns to be printed.  This
      also applies to vectors that are printed horizontally.
    }
    \item{printLevel}{0, the level of verbosity.  Larger values print
      more information.  Result depends on the optimizer.  Form
      \code{print.level} is also accepted by the methods for
      compatibility.}
    \item{storeParameters}{\code{FALSE}, whether to store and return the
      parameter
      values at each epoch.  If \code{TRUE}, the stored values
      can be retrieved with \code{\link{storedParameters}}-method.  The
      parameters are stored as a matrix with rows corresponding to the
      epochs and columns to the parameter components.
    }
    \item{storeValues}{\code{FALSE}, whether to store and return the objective
      function values at each epoch.  If \code{TRUE}, the stored values
      can be retrieved with \code{\link{storedValues}}-method.}
  }
}
\section{Methods}{
  \describe{
    \item{maxControl}{\code{(\dots)} creates a \dQuote{MaxControl} object.  The
      arguments must be in the form \code{option1 = value1, option2 =
	value2, ...}.  The options should be slot names, but the method
      also supports selected other parameter forms for compatibility reasons
      e.g. \dQuote{print.level} instead of \dQuote{printLevel}.
      In case there are more than one option with
      similar name, the last one overwrites the previous values.  This
      allows the user to override default parameters in the control
      list.  See example in \link{maxLik-package}.
    }
    \item{maxControl}{\code{(x = "MaxControl", \dots)} overwrites parameters
      of an existing \dQuote{MaxControl} object.  The \sQuote{\dots}
      argument must be in the form \code{option1 = value1, option2 =
	value2, ...}.  In case there are more than one option with
      similar name, only the last one is taken into account.  This
      allows the user to override default parameters in the control
      list. See example in \link{maxLik-package}.
    }
    \item{maxControl}{\code{(x = "maxim")} extracts \dQuote{MaxControl}
      structure from an estimated model}
    \item{show}{shows the parameter values}
  }
}
\section{Details}{
  Typically, the control options are supplied in the form of a list, in which
  case the corresponding default values are overwritten by the
  user-specified ones.  However, one may also create the control
  structure by \code{maxControl(opt1=value1, opt2=value2, ...)} and
  supply such value directly to the optimizer.  In this case the
  optimization routine takes all the values from the control object.
  }
\references{
  \itemize{
    \item Nelder, J. A. & Mead, R. A (1965) Simplex Method for Function
    Minimization \emph{The Computer Journal} \bold{7}, 308--313
    \item Marquardt, D. W. (1963) An Algorithm for Least-Squares Estimation of
    Nonlinear Parameters \emph{Journal of the Society for Industrial and
    Applied Mathematics} \bold{11}, 431--441 
  }
}
\author{
  Ott Toomet
}
\note{
  Several control parameters can also be supplied directly to the
  optimization routines.
}
\examples{
library(maxLik)
## Create a 'maxControl' object:
maxControl(tol=1e-4, sann_tmax=7, printLevel=2)

## Optimize quadratic form t(D) %*% W %*% D with p.d. weight matrix,
## s.t. constraints sum(D) = 1
quadForm <- function(D) {
   return(-t(D) \%*\% W \%*\% D)
}
eps <- 0.1
W <- diag(3) + matrix(runif(9), 3, 3)*eps
D <- rep(1/3, 3)
                        # initial values
## create control object and use it for optimization
co <- maxControl(printLevel=2, qac="marquardt", marquardt_lambda0=1)
res <- maxNR(quadForm, start=D, control=co)
print(summary(res))
## Now perform the same with no trace information
co <- maxControl(co, printLevel=0)
res <- maxNR(quadForm, start=D, control=co) # no tracing information
print(summary(res))  # should be the same as above
maxControl(res) # shows the control structure
}
\keyword{utilities}
