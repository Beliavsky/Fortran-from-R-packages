\name{pnvmix}
\alias{pnvmix}
\alias{pStudent}
\alias{pNorm}
\title{Distribution Function of Multivariate Normal Variance Mixtures}
\description{
  Evaluating multivariate normal variance mixture distribution functions
  (including Student \emph{t} and normal distributions).
}
\usage{
pnvmix(upper, lower = matrix(-Inf, nrow = n, ncol = d), qmix, rmix,
      loc = rep(0, d), scale = diag(d), standardized = FALSE,
      control = list(), verbose = TRUE, ...)

pStudent(upper, lower = matrix(-Inf, nrow = n, ncol = d), df, loc = rep(0, d),
         scale = diag(d), standardized = FALSE, control = list(), verbose = TRUE)
pNorm(upper, lower = matrix(-Inf, nrow = n, ncol = d), loc = rep(0, d),
      scale = diag(d), standardized = FALSE, control = list(), verbose = TRUE)
}
\arguments{
  \item{upper}{\eqn{(n, d)}-\code{\link{matrix}} of upper integration limits;
    each row represents a \eqn{d}-vector of upper integration limits.}
  \item{lower}{\eqn{(n, d)}-\code{\link{matrix}} of lower integration limits
    (componentwise less than or equal to \code{upper});
    each row represents a \eqn{d}-vector of lower integration limits.}
  \item{qmix, rmix}{specification of the mixing variable \eqn{W} via a quantile
  function (\code{qmix}) (recommended, see details below)
  *or* random number generator (\code{rmix}); see McNeil et
    al. (2015, Chapter 6) and Hintz et al. (2020). Supported are the following types of
    specification (see also the examples below):
    \describe{\item{\code{\link{character}}:}{\code{\link{character}} string
	specifying a supported distribution; currently available are
        \code{"constant"} (in which case \eqn{W = 1} and thus
	the multivariate normal distribution with mean vector
	\code{loc} and covariance matrix \code{scale} results),
	\code{"inverse.gamma"} (in which case \eqn{W} is
	inverse gamma distributed with shape and rate parameters
	\code{df}/2 and thus the multivariate
	Student \emph{t} distribution with \code{df} degrees of freedom
	(required to be provided via the ellipsis argument) results) and
	\code{"pareto"} (in which case \eqn{W} is
	Pareto distributed with scale equal to unity and shape equal to
	\code{alpha}, which needs to be provided via the ellipsis
        argument).}
      \item{\code{\link{list}}:}{\code{\link{list}} of length at least
	one, where the first component is a \code{\link{character}}
	string specifying the base name of a distribution whose
	quantile function or random number generator can be accessed via the prefix
	\code{"q"} and \code{"r"}, respectively.
	an example is \code{"exp"} for which \code{"qexp"} exists.
	If the list is of length larger than one, the remaining elements
	contain additional parameters of the distribution; for \code{"exp"},
        for example, this can be the parameter \code{rate}.}
      \item{\code{\link{function}}:}{\code{\link{function}}
	interpreted as the quantile function or random number generator of the mixing
	variable \eqn{W}.}
    }
  }
  \item{df}{positive degress of freedom; can also be \code{Inf} in which
    case the distribution is interpreted as the multivariate normal
    distribution with mean vector \code{loc} and covariance matrix
    \code{scale}.}
  \item{loc}{location vector of dimension \eqn{d}; this equals the mean
    vector of a random vector following the specified normal variance
    mixture distribution if and only if the latter exists.}
  \item{scale}{scale matrix (a covariance matrix entering the
    distribution as a parameter) of dimension \eqn{(d, d)};
    this equals the covariance matrix of a random vector following
    the specified normal variance mixture distribution divided by
    the expecation of the mixing variable \eqn{W} if and only if the
    former exists. \code{scale} is allowed to be singular in which case
    the distribution function of the singular normal variance mixture
    is returned.}
  \item{standardized}{\code{logical} indicating whether \code{scale} is
    assumed to be a correlation matrix.}
  \item{control}{\code{\link{list}} specifying algorithm specific
    parameters; see \code{\link{get_set_param}()}.}
  \item{verbose}{% *not* the same as in ?dnvmix()
    \code{\link{logical}}
    indicating whether a warning is thrown if the required precision
    \code{pnvmix.abstol} or \code{pnvmix.reltol} as specified in the
    \code{control} argument has not been reached; can also be an
    \code{\link{integer}} in which case 0 is \code{FALSE}, 1 is
    \code{TRUE} and 2 stands for producing a more verbose warning (for
    each set of provided integration bounds).}
  \item{\dots}{additional arguments (for example, parameters) passed
    to the underlying mixing distribution when \code{qmix} is a
    \code{\link{character}} string or \code{\link{function}}.}
}
\value{
  \code{pnvmix()}, \code{pStudent()} and \code{pNorm()} return a
  \code{\link{numeric}} \eqn{n}-vector with the computed probabilities
  and corresponding attributes \code{"abs. error"} and \code{rel. error}
  (error estimates of the RQMC estimator) and \code{"numiter"} (number of iterations).
}
\details{
  One should highlight that evaluating normal variance mixtures is a
  non-trivial tasks which, at the time of development of \pkg{nvmix},
  was not available in \R before, not even the special case of a
  multivariate Student \emph{t} distribution for non-integer degrees of
  freedom, which frequently appears in applications in finance,
  insurance and risk management after estimating such distributions.

  Note that the procedures call underlying C code. Currently, dimensions
  \eqn{d\ge 16510}{d >= 16510} are not supported for the default method
  \code{sobol}.

  Internally, an iterative randomized Quasi-Monte Carlo (RQMC) approach
  is used to estimate the probabilities. It is an iterative algorithm
  that evaluates the integrand at a point-set (with size as specified by
  \code{control$increment} in the \code{control} argument) in each
  iteration until the pre-specified absolute error tolerance
  \code{control$pnvmix.abstol} (or relative error tolerance
  \code{control$pnvmix.reltol} which is used only when
  \code{control$pnvmix.abstol = NA}) is reached.  The attribute
  \code{"numiter"} gives the number of such iterations needed.
  Algorithm specific parameters (such as the above mentioned
  \code{control$pnvmix.abstol}) can be passed as a list
  via \code{control}, see \code{\link{get_set_param}()} for more
  details. If specified error tolerances are not reached and
  \code{verbose = TRUE}, a warning is thrown.

  If provided \code{scale}
  is singular, \code{pnvmix()} estimates the correct probability but
  throws a warning if \code{verbose = TRUE}.

  It is recommended to supply a quantile function via \code{qmix}, if available,
  as in this case efficient RQMC methods are used to approximate the probability.
  If \code{rmix} is provided, internally used is
  plain MC integration, typically leading to slower convergence.
  If both \code{qmix} and \code{rmix} are provided, the latter
  is ignored.

  \code{pStudent()} and \code{pNorm()} are wrappers of
  \code{pnvmix(, qmix = "inverse.gamma", df = df)} and
  \code{pnvmix(, qmix = "constant")}, respectively.
  In the univariate case, the functions
  \code{\link{pt}()} and \code{\link{pnorm}()} are used.
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

  Genz, A. and Bretz, F. (1999).
  Numerical computation of multivariate t-probabilities with
  application to power calculation of multiple contrasts.
  \emph{Journal of Statistical Computation and Simulation} 63(4),
  103--117.

  Genz, A. and Bretz, F. (2002).
  Comparison of methods for the computation of multivariate \emph{t}
  probabilities.
  \emph{Journal of Computational and Graphical Statistics} 11(4),
  950--971.

  Genz, A. and Kwong, K. (2000).
  Numerical evaluation of singular multivariate normal distributions.
  \emph{Journal of Statistical Computation and Simulation} 68(1),
  1--21.
}
\seealso{
  \code{\link{dnvmix}()}, \code{\link{rnvmix}()}, \code{\link{fitnvmix}()},
  \code{\link{pgnvmix}()}, \code{\link{get_set_param}()}
}
\examples{
### Examples for pnvmix() ######################################################

## Generate a random correlation matrix in d dimensions
d <- 3
set.seed(157)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A \%*\% t(A))

## Evaluate a t_{1/2} distribution function
a <- -3 * runif(d) * sqrt(d) # random lower limit
b <-  3 * runif(d) * sqrt(d) # random upper limit
df <- 1.5 # note that this is *non-integer*
set.seed(123)
pt1 <- pnvmix(b, lower = a, qmix = "inverse.gamma", df = df, scale = P)

## Here is a version providing the quantile function of the mixing distribution
qmix <- function(u, df) 1 / qgamma(1-u, shape = df/2, rate = df/2)
mean.sqrt.mix <- sqrt(df) * gamma(df/2) / (sqrt(2) * gamma((df+1) / 2))
set.seed(123)
pt2 <- pnvmix(b, lower = a, qmix = qmix, df  = df, scale = P,
              control = list(mean.sqrt.mix = mean.sqrt.mix))

## Compare
stopifnot(all.equal(pt1, pt2, tol = 7e-4, check.attributes = FALSE))

## mean.sqrt.mix will be approximated by QMC internally if not provided,
## so the results will differ slightly.
set.seed(123)
pt3 <- pnvmix(b, lower = a, qmix = qmix, df = df, scale = P)
stopifnot(all.equal(pt3, pt1, tol = 7e-4, check.attributes = FALSE))

## Here is a version providing a RNG for the mixing distribution
## Note the significantly larger number of iterations in the attribute 'numiter'
## compared to when 'qmix' was provided (=> plain MC versus RQMC)
set.seed(123)
pt4 <- pnvmix(b, lower = a,
              rmix = function(n, df) 1/rgamma(n, shape = df/2, rate = df/2),
              df = df, scale = P)
stopifnot(all.equal(pt4, pt1, tol = 8e-4, check.attributes = FALSE))

## Case with missing data and a matrix of lower and upper bounds
a. <- matrix(rep(a, each = 4), ncol = d)
b. <- matrix(rep(b, each = 4), ncol = d)
a.[2,1] <- NA
b.[3,2] <- NA
pt <- pnvmix(b., lower = a., qmix = "inverse.gamma", df = df, scale = P)
stopifnot(is.na(pt) == c(FALSE, TRUE, TRUE, FALSE))

## Case where upper = (Inf,..,Inf) and lower = (-Inf,...,-Inf)
stopifnot(all.equal(pnvmix(upper = rep(Inf, d), qmix = "constant"), 1,
                    check.attributes = FALSE))

## An example with singular scale:
A <- matrix( c(1, 0, 0, 0,
               2, 1, 0, 0,
               3, 0, 0, 0,
               4, 1, 0, 1), ncol = 4, nrow = 4, byrow = TRUE)
scale <- A\%*\%t(A)
upper <- 2:5

pn <- pnvmix(upper, qmix = "constant", scale = scale) # multivariate normal
pt <- pnvmix(upper, qmix = "inverse.gamma", scale = scale, df = df) # multivariate t

stopifnot(all.equal(pn, 0.8581, tol = 1e-3, check.attributes = FALSE))
stopifnot(all.equal(pt, 0.7656, tol = 1e-3, check.attributes = FALSE))

## Evaluate a Exp(1)-mixture
## Specify the mixture distribution parameter
rate <- 1.9 # exponential rate parameter

## Method 1: Use R's qexp() function and provide a list as 'mix'
set.seed(42)
(p1 <- pnvmix(b, lower = a, qmix = list("exp", rate = rate), scale = P))

## Method 2: Define the quantile function manually (note that
##           we do not specify rate in the quantile function here,
##           but conveniently pass it via the ellipsis argument)
set.seed(42)
(p2 <- pnvmix(b, lower = a, qmix = function(u, lambda) -log(1-u)/lambda,
              scale = P, lambda = rate))

## Check
stopifnot(all.equal(p1, p2))


### Examples for pStudent() and pNorm() ########################################

## Evaluate a t_{3.5} distribution function
set.seed(271)
pt <- pStudent(b, lower = a, df = 3.5, scale = P)
stopifnot(all.equal(pt, 0.6180, tol = 7e-5, check.attributes = FALSE))

## Evaluate a normal distribution function
set.seed(271)
pn <- pNorm(b, lower = a, scale = P)
stopifnot(all.equal(pn, 0.7001, tol = 1e-4, check.attributes = FALSE))

## pStudent deals correctly with df = Inf:
set.seed(123)
p.St.dfInf <- pStudent(b, df = Inf, scale = P)
set.seed(123)
p.Norm <- pNorm(b, scale = P)
stopifnot(all.equal(p.St.dfInf, p.Norm, check.attributes = FALSE))
}
\keyword{distribution}