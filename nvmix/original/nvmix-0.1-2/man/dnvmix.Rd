\name{dnvmix}
\alias{dnvmix}
\alias{dStudent}
\alias{dNorm}
\title{Density of Multivariate Normal Variance Mixtures}
\description{
  Evaluating multivariate normal variance mixture densities
  (including Student \emph{t} and normal densities).
}
\usage{
dnvmix(x, qmix, loc = rep(0, d), scale = diag(d),
       factor = NULL, control = list(), log = FALSE, verbose = TRUE,...)

dStudent(x, df, loc = rep(0, d), scale = diag(d), factor = NULL,
         log = FALSE, verbose = TRUE, ...)
dNorm(x, loc = rep(0, d), scale = diag(d), factor = NULL,
      log = FALSE, verbose = TRUE, ...)
}
\arguments{
  \item{x}{\eqn{(n, d)}-\code{\link{matrix}} of evaluation points.}
  \item{qmix}{specification of the mixing variable \eqn{W}; see
  \code{\link{pnvmix}()} for details and examples.}
  \item{df}{positive degress of freedom; can also be \code{Inf} in which
    case the distribution is interpreted as the multivariate normal
    distribution with mean vector \code{loc} and covariance matrix
    \code{scale}).}
  \item{loc}{location vector of dimension \eqn{d}; this equals the mean
    vector of a random vector following the specified normal variance
    mixture distribution if and only if the latter exists.}
  \item{scale}{scale matrix (a covariance matrix entering the
    distribution as a parameter) of dimension \eqn{(d, d)};
    this equals the covariance matrix of a random vector following
    the specified normal variance mixture distribution divided by
    the expecation of the mixing variable \eqn{W} if and only if the
    former exists. Needs to be full rank for the density to exist.}
  \item{factor}{\eqn{(d, d)} lower triangular \code{\link{matrix}}
    such that \code{factor \%*\% t(factor)} equals \code{scale};
    note that for performance reasons, this property is not tested.
    If not provided, \code{factor} is internally determined via
    \code{t(\link{chol}())}.}
  \item{control}{\code{\link{list}} specifying algorithm specific
    parameters; see \code{\link{get_set_param}()}.}
  \item{log}{\code{\link{logical}} indicating whether the logarithmic
    density is to be computed.}
  \item{verbose}{\code{\link{logical}} indicating whether a
    warning is given if the required precision \code{abstol} has not
    been reached.}
  \item{\dots}{additional arguments (for example, parameters)
    passed to the underlying mixing distribution when \code{qmix}
    is a \code{\link{character}} string or \code{\link{function}}.}
}
\value{
  \code{dnvmix()}, \code{dStudent()} and \code{dNorm()} return a
  \code{\link{numeric}} \eqn{n}-vector with the computed (log-)density
  values and attributes \code{"abs. error"} and \code{"rel. error"}
  (containing the absolute and relative error
  estimates of the of the (log-)density) and \code{"numiter"}
  (containing the number of iterations).
}
\details{
  For the density to exist, \code{scale} must be full rank.
  Internally used is \code{factor}, so \code{scale} is not required
  to be provided if \code{factor} is given. The default factorization
  used to obtain \code{factor} is the Cholesky
  decomposition via \code{\link{chol}()}.

  \code{dStudent()} and \code{dNorm()} are wrappers of
  \code{dnvmix(, qmix = "inverse.gamma", df = df)} and
  \code{dnvmix(, qmix = "constant")}, respectively.
  In these cases, \code{dnvmix()} uses the analytical formulas for the
  density of a multivariate Student \emph{t} and normal distribution,
  respectively.

  Internally, an adaptive randomized Quasi-Monte Carlo (RQMC) approach
  is used to estimate the log-density. It is an iterative algorithm that
  evaluates the integrand at a randomized Sobol' point-set (default) in
  each iteration until the pre-specified error tolerance
  \code{control$dnvmix.reltol} in the \code{control} argument
  is reached for the log-density.  The attribute
  \code{"numiter"} gives the worst case number of such iterations needed
  (over all rows of \code{x}). Note that this function calls underlying
  C code.

  Algorithm specific parameters (such as above mentioned \code{control$dnvmix.reltol})
  can be passed as a \code{\link{list}} via the \code{control} argument,
  see \code{\link{get_set_param}()} for details and defaults.

  If the error tolerance cannot be achieved within \code{control$max.iter.rqmc}
  iterations and \code{fun.eval[2]} function evaluations, an additional
  warning is thrown if \code{verbose=TRUE}.
}
\author{Erik Hintz, Marius Hofert and Christiane Lemieux.}
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
}
\seealso{
  \code{\link{pnvmix}()}, \code{\link{rnvmix}()}, \code{\link{fitnvmix}()},
  \code{\link{get_set_param}()}.
}
\examples{
### Examples for dnvmix() ######################################################

## Generate a random correlation matrix in three dimensions
d <- 3
set.seed(271)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A \%*\% t(A))

## Evaluate a t_{3.5} density
df <- 3.5
x <- matrix(1:12/12, ncol = d) # evaluation points
dt1 <- dnvmix(x, qmix = "inverse.gamma", df = df, scale = P)
stopifnot(all.equal(dt1, c(0.013266542, 0.011967156, 0.010760575, 0.009648682),
                    tol = 1e-7, check.attributes = FALSE))

## Here is a version providing the quantile function of the mixing distribution
qW <- function(u, df) 1 / qgamma(1-u, shape = df/2, rate = df/2)
dt2 <- dnvmix(x, qmix = qW, df = df, scale = P)

## Compare
stopifnot(all.equal(dt1, dt2, tol = 5e-4, check.attributes = FALSE))

## Evaluate a normal density
dn <- dnvmix(x, qmix = "constant", scale = P)
stopifnot(all.equal(dn, c(0.013083858, 0.011141923, 0.009389987, 0.007831596),
                    tol = 1e-7, check.attributes = FALSE))

## Case with missing data
x. <- x
x.[3,2] <- NA
x.[4,3] <- NA
dt <- dnvmix(x., qmix = "inverse.gamma", df = df, scale = P)
stopifnot(is.na(dt) == rep(c(FALSE, TRUE), each = 2))

## Univariate case
x.. <- cbind(1:10/10) # (n = 10, 1)-matrix; vectors are considered rows in dnvmix()
dt1 <- dnvmix(x.., qmix = "inverse.gamma", df = df, factor = 1)
dt2 <- dt(as.vector(x..), df = df)
stopifnot(all.equal(dt1, dt2, check.attributes = FALSE))


### Examples for dStudent() and dNorm() ########################################

## Evaluate a t_{3.5} density
dt <- dStudent(x, df = df, scale = P)
stopifnot(all.equal(dt, c(0.013266542, 0.011967156, 0.010760575, 0.009648682),
                    tol = 1e-7, check.attributes = FALSE))

## Evaluate a normal density
x <- x[1,] # use just the first point this time
dn <- dNorm(x, scale = P)
stopifnot(all.equal(dn, 0.013083858, tol = 1e-7, check.attributes = FALSE))
}
\keyword{distribution}