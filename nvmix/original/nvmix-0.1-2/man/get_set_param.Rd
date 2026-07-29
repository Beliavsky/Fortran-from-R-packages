\name{get_set_param}
\alias{get_set_param}
\title{Algorithm-specific Parameters}
\description{
  Algorithm specific parameters for functionalities in the \code{nvmix}
  package, notably for \code{\link{fitnvmix}()}, \code{\link{dnvmix}()},
  \code{\link{pnvmix}()}, \code{\link{qnvmix}()}, \code{\link{pgammamix}()},
  \code{\link{dgammamix}()} and \code{\link{ES_nvmix}()} as well as the
  corresponding functions for grouped mixtures.
}
\usage{
get_set_param(control = list())
}
\arguments{
  \item{control}{\code{\link{list}} specifying algorithm specific
    parameters to beset; see below under details.}
}
\value{
  \code{get_set_param()} returns a \code{\link{list}} with more
  than 30 elements specifying algorithm specific parameters for the
  functions \code{\link{fitnvmix}()}, \code{\link{dnvmix}()},
  \code{\link{pnvmix}()}, \code{\link{qnvmix}()}, \code{\link{pgammamix}()},
  \code{\link{dgammamix}()} and \code{\link{ES_nvmix}()}, as well as the
  corresponding functions for grouped mixtures such as \code{\link{pgnvmix}()}
  and \code{\link{dgnvmix}()}.
  Parameter values passed to \code{get_set_param()} via the
  \code{control} argument overwrite the defaults; for parameters not
  specified in the \code{control} argument, the default values are being
  returned.
}
\details{
  For most functions in the \code{nvmix} package, internally, an
  iterative randomized Quasi-Monte Carlo (RQMC) approach is used to
  estimate probabilities, weights and (log-)densities. There are various
  parameters of underlying methods than can be changed.

  Algorithm specific parameters can be passed as a list via
  \code{control}.  It can contain any of the following:

  \describe{
    \item{For all algorithms:}{
      \describe{
        \item{\code{method}}{\code{\link{character}} string indicating the
          method to be used to compute the integral. Available are:
          \describe{
            \item{\code{"sobol"}:}{Sobol' sequence (default),}
            \item{\code{"ghalton"}:}{generalized Halton sequence,}
            \item{\code{"PRNG"}:}{plain Monte Carlo based on a pseudo-random
              number generator.}
          }
        }
        \item{\code{increment}}{\code{\link{character}} string indicating how
          the sample size should be increased in each iteration. Available are:
          \describe{
            \item{\code{"doubling"}:}{next iteration has as many sample points
              as all the previous iterations combined,}
            \item{\code{"num.init"}:}{all iterations use an additional
              \code{fun.eval[1]}-many points (default for most functions).}
          }
        }
        \item{\code{CI.factor}}{multiplier of the Monte Carlo confidence interval
          bounds. The algorithm runs until \code{CI.factor} times the estimated
          standard error is less than \code{abstol} or \code{reltol} (whichever
          is provided). If \code{CI.factor = 3.5} (the default), one can expect
          the actual absolute error to be less than \code{abstol} in
	  99.9\% of the cases.}
        \item{\code{fun.eval}}{\code{\link{numeric}(2)} providing the size of
          the first point set to be used to estimate integrals
	  (typically a power of 2) and the maximal number of function
          evaluations. \code{fun.eval} defaults to \code{c(2^7, 1e12)}.}
        \item{\code{max.iter.rqmc}}{\code{\link{numeric}}, providing the maximum
          number of iterations allowed in the RQMC approach; the default is \code{15}
          if \code{increment = "doubling"} and \code{1250} otherwise.}
        \item{\code{B}}{number of randomizations for obtaining an error estimate
          in the RQMC approach; the default is \code{15}.}
      }
    }
    \item{For \code{pnvmix()} and \code{pgnvmix()}:}{
      \describe{
	\item{\code{pnvmix.abstol}, \code{pnvmix.reltol}}{non-negative \code{\link{numeric}}
          providing the relative/absolute precision required for the distribution
          function. Relative precision via \code{pnvmix.reltol} is only used
          when \code{pnvmix.abstol = NA}; in all other cases, absolute precision
          will be used. \code{pnvmix.abstol} defaults to \code{1e-3}.
          If \code{pnvmix.abstol = 0} and \code{pnvmix.reltol = 0}, the algorithm
          will typically run until the total number of function evaluations
          exceeds \code{fun.eval[2]} or until the total number of iterations exeeds
          \code{max.iter.rqmc}, whichever happens first.
          If \eqn{n > 1} (so \code{upper} has more than
          one row), the algorithm runs until the precision requirement is reached
          for all \eqn{n} probability estimates.}
	\item{\code{mean.sqrt.mix}}{expectation of the square root \eqn{\sqrt(W)}
          of the mixing variable \eqn{W}. If \code{NULL}, it will be estimated via
          QMC; this is only needed for determining the reordering of the
          integration bounds, so a rather crude approximation is fine.}
	\item{\code{precond}}{\code{\link{logical}} indicating whether
          preconditioning is applied, that is, reordering of the integration
          variables. If \code{TRUE}, integration limits \code{lower}, \code{upper}
          as well as \code{scale}
          are internally re-ordered in a way such that the overall variance of the
          integrand is usually smaller than with the original ordering; this
          usually leads smaller run-times.}
	\item{\code{cholesky.tol}}{non-negative numeric providing lower threshold
          for non-zero elements in the computation of the cholesky factor: If
          calculated \eqn{C(i,i)^2 < | cholesky.tol * Scale(i,i)|}, the diagonal
          element (and all other elements in column \eqn{i}) of the cholesky factor
          \eqn{C} are set to zero, yielding a singular matrix. \code{cholesky.tol}
          defaults to \code{1e-9}.}
      }
    }
    \item{For \code{dnvmix()} and \code{dgnvmix()}:}{
      \describe{
	\item{\code{dnvmix.reltol}, \code{dnvmix.abstol}}{non-negative
	\code{\link{numeric}} providing the relative/absolute precision for the *log-*
	  density required. Absolute precision via \code{dnvmix.abstol}
	  is only used when \code{dnvmix.reltol = NA}; in all other
	  cases, relative precision will be used. \code{dnvmix.reltol}
	  defaults to \code{1e-2}.
	  If \code{dnvmix.reltol=0} and \code{dnvmix.abstol=0}, the algorithm
	  will typically run until the total number of function evaluations exceeds
	  \code{fun.eval[2]} or until the total number of iterations exeeds
	  \code{max.iter.rqmc}, whichever happens first.
	  If \eqn{n > 1} (so \code{x} has more than one row), the algorithm runs until
	  the precision requirement is reached for all \eqn{n} log-density estimates.}
	\item{\code{dnvmix.doAdapt}}{\code{\link{logical}} indicating if an adaptive
	integration procedure shall be used that only samples in relevant subdomains
	of the mixing variable; defaults to \code{TRUE}.}
	\item{\code{dnvmix.max.iter.rqmc.pilot}}{\code{\link{numeric}}, providing
          the maximum number of unstratified, non-adaptive pilot runs the internal
          integration procedure performs. Defaults to \code{6}.}
	\item{\code{dnvmix.tol.int.lower}, \code{dnvmix.order.lower}}{
          both \code{\link{numeric}} and nonnegative. RQMC integration
          is only performed where the integrand is > than the maximum of
          \code{dnvmix.tol.int.lower} and \eqn{10^{-c} g_{max}}, where
          \eqn{g_{max}} is the theoretical maximum of the integrand and \eqn{c}
          is the specified \code{dnvmix.order.lower}. Default to \code{1e-100}
	  and \code{5}, respectively.}
	\item{\code{dnvmix.tol.bisec}}{\code{\link{numeric}} \code{vector} of
          length 3 specifying bisection tolerances in the adaptive RQMC algorithm.
          First/second/third element specify the tolerance on \eqn{u}, \eqn{W} and
          the log-integrand and default to \code{1e-6}, \code{1e-1} and
          \code{1e-1}, respectively.}
	\item{\code{dnvmix.max.iter.bisec}}{\code{\link{numeric}}, maximum number
          of iterations in the internal bisection procedure to find good cutting
          points allowed, defaults to \code{15}.}
	\item{\code{dnvmix.tol.stratlength}}{\code{\link{numeric}}, nonnegative.
          If the stratum found by the adaptive integration method has length >
          \code{dnvmix.tol.stratlength} RQMC integration is used there; otherwise
          a crude approximation. Defaults to \code{1e-50}.}
      }
    }
    \item{For \code{fitnvmix()}:}{
      \describe{
        \item{\code{ECMEstep}}{\code{\link{logical}}, if \code{TRUE} (default),
          ECME iteration is performed;
          if \code{FALSE}, no ECME step is performed so that \code{fitnvmix()}
          performs between zero and two optimizations over \eqn{nu}, depending
          on \code{laststep.do.nu} and whether \code{nu.init} was provided.}
        \item{\code{ECMEstep.do.nu}}{\code{\link{logical}}, if \code{TRUE}
          (default), the likelihood is maximized over \eqn{nu} in each ECME
          iteration; if \code{FALSE}, this step is omitted.}
        \item{\code{laststep.do.nu}}{\code{\link{logical}}, if \code{TRUE}
          another last maximization of the likelihood over \eqn{nu} is
          performed using all  observations after the ECME iterations.
          Only makes sense if either \code{ECMEstep.do.nu=FALSE} or
          if \code{size.subsample} is smaller than the number of observations.
          Defaults to \code{FALSE}.}
        \item{\code{resample}}{\code{\link{logical}}, if \code{TRUE}, a different
          subsample of \code{x} is taken in each optimization over \eqn{nu} in
          the ECME iterations. Only relevant when \code{size.subsample} is
          smaller than the number of observations. Defaults to \code{FALSE}.}
        \item{\code{ECME.miniter}, \code{ECME.maxiter}}{\code{\link{numeric}}
        positive, minimum and maximum number
          of ECME iterations. Default to \code{5} and \code{200}, respectively.}
        \item{\code{max.iter.locscaleupdate}}{\code{\link{numeric}} positive.
          Maximum number of location-scale updates (while helding \eqn{nu} fixed)
          in each individual ECME iteration; defaults to \code{50}.}
        \item{\code{weights.reltol}}{\code{numeric} non-negative. Relative
          tolerance to estimate internal weights used to update
	  \eqn{loc} and \eqn{scale} estimates in the ECME iterations. Defaults to
	  \code{1e-2}.}
        \item{\code{weights.interpol.reltol}}{\code{numeric}
	  non-negative. Some weights can be obtained by interpolating
	  previously calculated weights; if the maximal relative
	  interpolation error is smaller than
	  \code{weights.interpol.reltol}, this is done. Defaults to \code{1e-2}.}
        \item{\code{ECME.rel.conv.tol}}{\code{\link{numeric}(3)} vector
          specifying relative convergence tolerances for \code{loc}, \code{scale}
          and \code{nu} (in this order). Defaults to \code{c(1e-2, 1e-2, 1e-3)}.}
        \item{\code{control.optim}}{\code{\link{list}} of control parameters
          passed to the underlying \code{optim} in the initial step as well as
          in the ECME iterations. See \code{\link{optim}()} for details; defaults to
          \code{list(maxit=75)}.}
        \item{\code{control.optim.laststep}}{like \code{control.optim};
          this list of control arguments is passed to \code{optim} in the last-step.
          Only relevant when \code{laststep.do.nu = TRUE} and defaults to
          \code{list()} (so no defaults of \code{\link{optim}()} changed).}
      }
    }
    \item{For \code{qnvmix()}:}{
      \describe{
        \item{\code{max.iter.newton}}{\code{\link{numeric}}, maximum number
        of Newton iterations allowed to approximate the quantile; defaults
        to \code{45}.}
        \item{\code{newton.conv.abstol}}{\code{\link{numeric}}, convergence
        tolerance for the Newton proceudre; convergence is detected once
        the difference of estimated quantiles in two successive iterations
        is smaller than \code{newton.conv.abstol}; defaults to \code{5e-4}.}
        \item{\code{newton.df.reltol}}{\code{\link{numeric}}, relative error
        tolerance for estimating the univariate distribution function;
        defaults to \code{2.5e-4}.}
        \item{\code{newton.logdens.abstol}}{\code{\link{numeric}}, absolute error
        tolerance for the internal estimation of the log-density needed for the
        update; defaults to \code{1e-2}.}
        \item{\code{newton.df.max.iter.rqmc}}{\code{\link{numeric}}, maximum
        number of iterations to estimate the univariate distribution function
        required in the Newton update; defaults to \code{350}. Note that
        internally used is \code{increment = "doubling"}, no matter what.}
       }
    }
    \item{For \code{qqplot_maha()}:}{
      \describe{
        \item{\code{qqplot.df.reltol}}{\code{\link{numeric}}, with the same
        meaning as \code{newton.df.reltol} for the function \code{qnvmix()}.
        Defaults to \code{5e-3}.}
       }
    }
    \item{For \code{ES_nvmix()}:}{
      \describe{
        \item{\code{riskmeasures.abstol, riskmeasures.reltol}}{
        \code{\link{numeric}}, absolute or relative error tolerance for
        estimating riskmeasures, notably for \code{ES_nvmix()}. By default,
         \code{riskmeasures.reltol=5e-2} and \code{riskmeasures.abstol=NA},
         so that a relative tolerance is used.}
       }
    }
  }
  Care should be taken when changing algorithm specific parameters, notably
  tolerances, as the accuracy of the result is heavily influenced by those.
}
\author{Erik Hintz, Marius Hofert and Christiane Lemieux}
\references{
  Hintz, E., Hofert, M. and Lemieux, C. (2020),
  Grouped Normal Variance Mixtures.
  \emph{Risks} 8(4), 103.

  Hintz, E., Hofert, M. and Lemieux, C. (2021),
  Normal variance mixtures: Distribution, density and parameter estimation.
  \emph{Computational Statistics and Data Analysis} 157C, 107175.

  Hintz, E., Hofert, M. and Lemieux, C. (2022),
  Multivariate Normal Variance Mixtures in \R: The \R Package nvmix.
  \emph{Journal of Statistical Software}, \doi{10.18637/jss.v102.i02}.
}
\seealso{
  \code{\link{fitnvmix}()}, \code{\link{dnvmix}()},
  \code{\link{pnvmix}()}, \code{\link{qnvmix}()}, \code{\link{pgammamix}()},
  \code{\link{dgammamix}()}, \code{\link{ES_nvmix}()}
}
\examples{
get_set_param() # obtain defaults
}
\keyword{distribution}