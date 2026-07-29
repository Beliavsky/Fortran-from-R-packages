\name{fitnvmix}
\alias{fitnvmix}
\alias{fitStudent}
\alias{fitNorm}
\title{Fitting Multivariate Normal Variance Mixtures}
\description{
  Functionalities for fitting multivariate normal variance mixtures (in
  particular also Multivariate \emph{t} distributions) via an ECME algorithm.
}
\usage{
fitnvmix(x, qmix, mix.param.bounds, nu.init = NA, loc = NULL, scale = NULL,
         init.size.subsample = min(n, 100), size.subsample = n,
         control = list(), verbose = TRUE)

fitStudent(x, loc = NULL, scale = NULL, mix.param.bounds = c(1e-3, 1e2), ...)
fitNorm(x)
}
\arguments{
  \item{x}{\eqn{(n, d)}-data \code{\link{matrix}}.}
  \item{qmix}{specification of the mixing variable \eqn{W}; see McNeil et
    al. (2015, Chapter 6). Supported are the following types of
    specification (see also the examples below):
    \describe{
      \item{\code{\link{character}}:}{\code{\link{character}} string
	specifying a supported distribution; currently available are
        \code{"constant"} (in which case \eqn{W = 1} and thus
	the multivariate normal distribution with mean vector
	\code{loc} and covariance matrix \code{scale} results),
	\code{"inverse.gamma"} (in which case \eqn{W} is
	inverse gamma distributed with shape and rate parameters
	\code{df}/2 and thus the multivariate
	Student \emph{t} distribution with \code{df} degrees of freedom
	results) and \code{"pareto"} (in which case \eqn{W} is
	Pareto distributed with scale equal to unity and shape equal to
	\code{alpha}).
      }
      \item{\code{\link{function}}:}{\code{\link{function}}
	interpreted as the quantile function of the mixing
	variable \eqn{W}. In this case, \code{qmix} \emph{must} have the
	form \code{qmix = function(u, nu)}, where the argument \code{nu}
	corresponds to the parameter (vector) specifying the distribution
	of the mixing variable.
      }
    }
  }
  \item{mix.param.bounds}{either \code{\link{numeric}(2)} or a
    \code{\link{matrix}} with two columns. The first/second column
    corresponds to the lower/upper bound of \eqn{nu_i}, the ith
    component of the parameter vector \eqn{nu} of the mixing variable
    \eqn{W}. All elements need to be finite, numeric values.
    Note: The algorithm tends to converge quicker if the parameter
    ranges supplied are smaller.}
  \item{nu.init}{either \code{NA} or an initial guess for the parameter
    (vector) \eqn{nu}. In the former case an initial estimate is
    calculated by the algorithm. If \code{nu.init} is provided, the
    algorithm often converges faster; the better the starting value,
    the faster convergence.}
  \item{loc}{\eqn{d}-\code{\link{vector}}; if provided, taken as the 'true'
    location vector in which case \code{loc} is not estimated.}
  \item{scale}{positive definite \eqn{(d, d)}-\code{\link{matrix}}; if provided,
    taken as the 'true' scale matrix in which case \code{scale} is not estimated.}
  \item{init.size.subsample}{\code{\link{numeric}}, non-negative,
    giving the sub-samplesize used to get an initial estimate for
    \eqn{nu}. Only used if \code{is.na(nu.init)}, otherwise ignored.}
  \item{size.subsample}{\code{\link{numeric}}, non-negative, specifying
    the size of the subsample that is being used in the ECME iterations
    to optimize the log-likelihood over \eqn{nu}. Defaults to \code{n},
    so that the full sample is being used. Decreasing this number can
    lead to faster run-times (as fewer log-densities need to be
    estimated) but also to an increase in bias and variance.}
  \item{control}{\code{\link{list}} specifying algorithm specific
    parameters; see below under 'Details' and \code{\link{get_set_param}()}.}
  \item{verbose}{\code{\link{numeric}} or \code{\link{logical}} (in
    which case it is converted to \code{numeric}) specifying the
    amount of tracing to be done. If \code{0} or \code{FALSE}, neither tracing
    nor warnigns are communicated; if \code{1}, only warnigns are communicated,
    if \code{2} or \code{3}, warnings and (shorter or longer) tracing information is
    provided.}
  \item{\dots}{additional arguments passed to the underlying
    \code{fitnvmix()}.}
}
\value{
  The function \code{fitnvmix()} returns an S3 object of
  \code{\link{class}} \code{"fitnvmix"}, basically a \code{\link{list}}
  which contains, among others, the components
  \describe{
    \item{\code{nu}}{Estimated mixing parameter (vector) \code{nu}.}
    \item{\code{loc}}{Estimated or provided \code{loc} vector.}
    \item{\code{scale}}{Estimated or provided \code{scale} matrix.}
    \item{\code{max.ll}}{Estimated log-likelihood at reported estimates.}
    \item{\code{x}}{Input data matrix \code{x}.}
   }
   The methods \code{print()}, \code{summary()} and \code{plot()} are defined
   for the class \code{"fitnvmix"}.

  \code{fitStudent()} is a wrapper to \code{fitnvmix()} for parameter
  estimation of multivariate Student \emph{t} distributions; it also
  returns an S3 object of \code{\link{class}} \code{"fitnvmix"} where
  the fitted degrees of freedom are called \code{"df"} instead of
  \code{"nu"} (to be consistent
  with the other wrappers for the Student \emph{t} distributions).

  \code{fitNorm()} just returns a \code{\link{list}} with components
  \code{loc} (columnwise sample means) and \code{scale} (sample
  covariance matrix).
}
\details{
  The function \code{fitnvmix()} uses an ECME algorithm to approximate the
  MLEs of the parameters \code{nu}, \code{loc} and \code{scale} of a
  normal variance mixture specified by \code{qmix}. The underlying
  procedure successively estimates \code{nu} (with given \code{loc} and
  \code{scale}) by maximizing the likelihood which is approximated by
  \code{\link{dnvmix}()} (unless \code{qmix} is a \code{character}
  string, in which case analytical formulas for the log-densities are
  used) and \code{scale} and \code{loc} (given \code{nu}) using weights
  (which again need to be approximated) related to the posterior
  distribution, details can be found in the first reference below.

  It should be highlighted that (unless unless \code{qmix} is a
  \code{character} string), every log-likelihood and every weight needed
  in the estimation is numerically approximated via RQMC methods. For
  large dimensions and sample sizes this procedure can therefore be
  slow.

  Various tolerances and convergence criteria can be changed by the user
  via the \code{control} argument. For more details, see
  \code{\link{get_set_param}()}.
}
\author{Erik Hintz, Marius Hofert and Christiane Lemieux}
\references{
  Hintz, E., Hofert, M. and Lemieux, C. (2021),
  Normal variance mixtures: Distribution, density and parameter estimation.
  \emph{Computational Statistics and Data Analysis} 157C, 107175.

  Hintz, E., Hofert, M. and Lemieux, C. (2022),
  Multivariate Normal Variance Mixtures in \R: The \R Package nvmix.
  \emph{Journal of Statistical Software}, \doi{10.18637/jss.v102.i02}.

  McNeil, A. J., Frey, R. and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.

  Liu, C. and Rubin, D. (1994).
  The ECME algorithm: a simple extension of EM and ECM with faster monotone
  convergence.
  \emph{Biometrika} 81(4),
  633--648.
}
\seealso{
  \code{\link{dnvmix}()}, \code{\link{rnvmix}()}, \code{\link{pnvmix}()},
  \code{\link{qqplot_maha}()}, \code{\link{get_set_param}()}.
}
\examples{
## Sampling parameters
set.seed(274) # for reproducibility
nu               <- 2.8 # parameter used to sample data
d                <- 4 # dimension
n                <- 75 # small sample size to have examples run fast
loc              <- rep(0, d) # location vector
A                <- matrix(runif(d * d), ncol = d)
diag_vars        <- diag(runif(d, min = 2, max = 5))
scale            <- diag_vars \%*\% cov2cor(A \%*\% t(A)) \%*\% diag_vars # scale matrix
mix.param.bounds <- c(1, 5) # nu in [1, 5]

### Example 1: Fitting a multivariate t distribution ###########################

if(FALSE){
    ## Define 'qmix' as the quantile function of an IG(nu/2, nu/2) distribution
    qmix <- function(u, nu) 1 / qgamma(1 - u, shape = nu/2, rate = nu/2)
    ## Sample data using 'rnvmix'
    x <- rnvmix(n, qmix = qmix, nu = nu, loc = loc, scale = scale)
    ## Call 'fitvnmix' with 'qmix' as a function (so all densities/weights are estimated)
    (MyFit11 <- fitnvmix(x, qmix = qmix, mix.param.bounds = mix.param.bounds))
    ## Call 'fitnvmix' with 'qmix = "inverse.gamma"' in which case analytical formulas
    ## for weights and densities are used:
    (MyFit12 <- fitnvmix(x, qmix = "inverse.gamma",
                         mix.param.bounds = mix.param.bounds))
    ## Alternatively, use the wrapper 'fitStudent()'
    (MyFit13 <- fitStudent(x))
    ## Check
    stopifnot(all.equal(MyFit11$nu, MyFit12$nu, tol = 5e-2),
              all.equal(MyFit11$nu, MyFit13$nu, tol = 5e-2))
    ## Can also provide 'loc' and 'scale' in which case only 'nu' is estimated
    (MyFit13 <- fitnvmix(x, qmix = "inverse.gamma", mix.param.bounds = mix.param.bounds,
                         loc = loc, scale = scale))
    (MyFit14 <- fitStudent(x, loc = loc, scale = scale))
    stopifnot(all.equal(MyFit13$nu, MyFit14$df, tol = 1e-6))
}

### Example 2: Fitting a Pareto mixture ########################################

## Define 'qmix' as the quantile function of a Par(nu, 1) distribution
qmix <- function(u, nu) (1-u)^(-1/nu)
## Sample data using 'rnvmix':
x <- rnvmix(n, qmix = qmix, nu = nu, loc = loc, scale = scale)
## Call 'fitvnmix' with 'qmix' as function (=> densities/weights estimated)
(MyFit21 <- fitnvmix(x, qmix = qmix, mix.param.bounds = mix.param.bounds))
## Call 'fitnvmix' with 'qmix = "pareto"' in which case an analytical formula
## for the density is used
(MyFit22 <- fitnvmix(x, qmix = "pareto", mix.param.bounds = mix.param.bounds))
stopifnot(all.equal(MyFit21$nu, MyFit22$nu, tol = 5e-2))

}
\keyword{distribution}