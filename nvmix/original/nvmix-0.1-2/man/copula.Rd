\name{copula}
\alias{dnvmixcopula}
\alias{pnvmixcopula}
\alias{rnvmixcopula}
\alias{dStudentcopula}
\alias{pStudentcopula}
\alias{rStudentcopula}
\alias{pgStudentcopula}
\alias{dgStudentcopula}
\alias{rgStudentcopula}
\alias{fitgStudentcopula}
\alias{fitStudentcopula}
\title{Functionalities for Normal Variance Mixture Copulas}
\description{
  Evaluate the density / distribution function of normal variance mixture
  copulas (including Student \emph{t} and normal copula) and generate vectors
  of random variates from normal variance mixture copulas.
}
\usage{
dnvmixcopula(u, qmix, scale = diag(d), factor = NULL, control = list(),
             verbose = FALSE, log = FALSE, ...)
pnvmixcopula(upper, lower = matrix(0, nrow = n, ncol = d), qmix, scale = diag(d),
             control = list(), verbose = FALSE, ...)
rnvmixcopula(n, qmix, scale = diag(2), factor = NULL,
             method = c("PRNG", "sobol", "ghalton"), skip = 0,
             control = list(), verbose = FALSE, ...)

dStudentcopula(u, df, scale = diag(d), factor = NULL, log = FALSE, verbose = TRUE)
pStudentcopula(upper, lower = matrix(0, nrow = n, ncol = d), df, scale = diag(d),
               control = list(), verbose = TRUE)
rStudentcopula(n, df, scale = diag(2), method = c("PRNG", "sobol", "ghalton"),
               skip = 0)

pgStudentcopula(upper, lower = matrix(0, nrow = n, ncol = d), groupings = 1:d,
                df, scale = diag(d), control = list(), verbose = TRUE)
dgStudentcopula(u, groupings = 1:d, df, scale = diag(d), factor = NULL,
                factor.inv = NULL, control = list(), verbose = TRUE, log = FALSE)
rgStudentcopula(n, groupings = 1:d, df, scale = diag(2), factor = NULL,
                method = c("PRNG", "sobol", "ghalton"), skip = 0)

fitgStudentcopula(x, u, df.init = NULL, scale = NULL, groupings = rep(1, d),
                  df.bounds = c(0.5, 30), fit.method = c("joint-MLE",
                  "groupewise-MLE"), control = list(), verbose = TRUE)
fitStudentcopula(u, fit.method = c("Moment-MLE", "EM-MLE", "Full-MLE"),
                 df.init = NULL, df.bounds = c(0.1, 30), control = list(),
                 verbose = TRUE)

}
\arguments{
  \item{u}{\eqn{(n, d)}-\code{\link{matrix}} of evaluation points or data;
    Have to be in \eqn{(0,1)}.}
  \item{upper, lower}{\eqn{(n, d)}-\code{\link{matrix}} of upper/lower evaluation
  limits. Have to be in \eqn{(0,1)}.}
  \item{n}{sample size \eqn{n} (positive integer).}
  \item{qmix}{specification of the mixing variable \eqn{W}; see
  \code{\link{pnvmix}()} for the ungrouped and \code{\link{pgnvmix}()} for the
  grouped case.}
  \item{groupings}{see \code{\link{pgnvmix}()}.}
  \item{df}{positive degress of freedom; can also be \code{Inf} in which
    case the copula is interpreted as the Gaussian copula.}
  \item{scale}{scale matrix (a covariance matrix entering the
    distribution as a parameter) of dimension \eqn{(d, d)} (defaults to
    \eqn{d = 2});
    this equals the covariance matrix of a random vector following
    the specified normal variance mixture distribution divided by
    the expecation of the mixing variable \eqn{W} if and only if the
    former exists.
    Note that \code{scale} must be positive definite; sampling from
    singular ungrouped normal variance mixtures can be achieved by providing
    \code{factor}.}
  \item{factor}{% more general here than for dnvmix() or pnvmix()
    \eqn{(d, k)}-\code{\link{matrix}} such that
    \code{factor \%*\% t(factor)} equals \code{scale}; the non-square
    case \eqn{k \neq d}{k != d} can be used to sample from singular
    normal variance mixtures. For \code{dnvmixcopula()}, this has to be
    a square matrix. Note that this notation coincides with
    McNeil et al. (2015, Chapter 6). If not provided, \code{factor} is
    internally determined via \code{\link{chol}()} (and multiplied
    from the right to an \eqn{(n, k)}-matrix of independent
    standard normals to obtain a sample from a multivariate normal
    with zero mean vector and covariance matrix \code{scale}).}
  \item{factor.inv}{inverse of \code{factor}; if not provided, computed via
    \code{solve(factor)}.}
  \item{method}{see \code{\link{rnvmix}()}.}
  \item{skip}{see \code{\link{rnvmix}()}.}
  \item{df.init}{\code{\link{NULL}} or vector with initial estimates for \code{df};
  can contain NAs.}
  \item{df.bounds}{\eqn{2}-\code{\link{vector}} with the lower/upper bounds on
  the degree-of-freedom parameter for the fitting.}
  \item{fit.method}{\code{\link{character}} indicating which fitting method is to
   be used; see details below.}
  \item{x}{\eqn{(n, d)}-\code{\link{matrix}} data matrix of which the underlying
  copula is to be estimated. See also details below.}
  \item{control}{\code{\link{list}} specifying algorithm specific
    parameters; see \code{\link{get_set_param}()}.}
  \item{verbose}{
    \code{\link{logical}}
    indicating whether a warning is given if the required precision
    \code{abstol} has not been reached.}
  \item{log}{\code{\link{logical}} indicating whether the logarithmic
    density is to be computed.}
  \item{\dots}{additional arguments (for example, parameters) passed to the
    underlying mixing distribution when \code{rmix} or \code{qmix} is a
    \code{\link{character}} string or \code{\link{function}}.}
}
\value{
  The values returned by \code{dnvmixcopula()}, \code{rnvmixcopula()} and
  \code{pnvmixcopula()} are similar to the ones returned by their
  non-copula alternatives \code{\link{dnvmix}()}, \code{\link{rnvmix}()}
  and \code{\link{pnvmix}()}.

  The function \code{fitgStudentcopula()} returns an S3 object of
  \code{\link{class}} \code{"fitgStudentcopula"}, basically a \code{\link{list}}
  which contains, among others, the components
  \describe{
    \item{\code{df}}{Estimated degrees-of-freedom for each group.}
    \item{\code{scale}}{Estimated or provided \code{scale} matrix.}
    \item{\code{max.ll}}{Estimated log-likelihood at reported estimates.}
    \item{\code{df.init}}{Initial estimate for the degrees-of-freedom.}
   }
   The methods \code{print()} and \code{summary()} are defined for the class
   \code{"fitgStudentcopula"}.
}
\details{
  Functionalities for normal variance mixture copulas provided here
  essentially call \code{\link{pnvmix}()}, \code{\link{dnvmix}()} and
  \code{\link{rnvmix}()} as well as \code{\link{qnvmix}()}, see their
  documentations for more details.

  We remark that computing normal variance mixtures is a challenging
  task; evaluating normal variance mixture copulas additionally requires
  the approximation of a univariate quantile function so that for large
  dimensions and sample sizes, these procedures can be fairly slow. As
  there are approximations on many levels, reported error estimates for
  the copula versions of \code{pnvmix}() and \code{dnvmix}() can be
  flawed.

  The functions \code{[d/p/r]Studentcopula()} are user-friendly wrappers for
  \code{[d/p/r]nvmixcopula(, qmix = "inverse.gamma")}, designed for the imporant
  case of a t copula with degrees-of-freedom \code{df}.

  The function \code{fitgStudentcopula()} can be used to estimate the matrix
  \code{scale} and the degrees-of-freedom for grouped t-copulas. The matrix
  \code{scale}, if not provided, is estimated non-parametrically. Initial values
  for the degrees-of-freedom are estimated for each group separately (by fitting
  the corresponding marginal t copula). Using these initial values, the joint
  likelihood over all (\code{length(unique(groupings))}-many) degrees-of-freedom
  parameters is optimized via \code{\link{optim}()}. For small dimensions,
  the results are satisfactory but the optimization becomes extremely challenging
  when the dimension is large, so care should be taking when interpreting the
  results.
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

  McNeil, A. J., Frey, R. and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.

  Luo, X. and Shevchenko, P. (2010).
  The \emph{t} copula with multiple parameters of degrees of freedom: bivariate
  characteristics and application to risk management.
  \emph{Quantitative Finance} 10(9),
  1039-1054.

  Daul, S., De Giorgi, E. G., Lindskog, F. and McNeil, A (2003).
  The grouped \emph{t} copula with an application to credit risk.
  \emph{Available at SSRN 1358956}.
}
\seealso{
  \code{\link{dnvmix}()}, \code{\link{pnvmix}()}, \code{\link{qnvmix}()},
  \code{\link{rnvmix}()}
}
\examples{
## Generate a random correlation matrix in d dimensions
d <- 2 # dimension
set.seed(42) # for reproducibility
rho <- runif(1, min = -1, max = 1)
P <- matrix(rho, nrow = d, ncol = d) # build the correlation matrix P
diag(P) <- 1
## Generate two random evaluation points:
u <- matrix(runif(2*d), ncol = d)
## We illustrate using a t-copula
df = 2.1
## Define quantile function which is inverse-gamma here:
qmix. <- function(u) 1/qgamma(1-u, shape = df/2, rate = df/2)


### Example for dnvmixcopula() ####################################################

## If qmix = "inverse.gamma", dnvmix() calls qt and dt:
d1 <- dnvmixcopula(u, qmix = "inverse.gamma", scale = P, df = df)
## Same can be obtained using 'dStudentcopula()'
d2 <- dStudentcopula(u, scale = P, df = df)
stopifnot(all.equal(d1, d2))
## Use qmix. to force the algorithm to use a rqmc procedure:
d3 <- dnvmixcopula(u, qmix = qmix., scale = P)
stopifnot(all.equal(d1, d3, tol = 1e-3, check.attributes = FALSE))


### Example for pnvmixcopula() ####################################################

## Same logic as above:
p1 <- pnvmixcopula(u, qmix = "inverse.gamma", scale = P, df = df)
p2 <- pnvmixcopula(u, qmix = qmix., scale = P)
stopifnot(all.equal(p1, p2, tol = 1e-3, check.attributes = FALSE))


### Examples for rnvmixcopula() ###################################################

## Draw random variates and compare
n <- 60
set.seed(1)
X  <- rnvmixcopula(n, qmix = "inverse.gamma", df = df, scale = P) # with scale
set.seed(1)
X. <- rnvmixcopula(n, qmix = "inverse.gamma", df = df, factor = t(chol(P))) # with factor
stopifnot(all.equal(X, X.))


### Example for the grouped case ##################################################

d <- 4 # dimension
set.seed(42) # for reproducibility
P <- matrix(runif(1, min = -1, max = 1), nrow = d, ncol = d) # build the correlation matrix P
diag(P) <- 1
groupings <- c(1, 1, 2, 2) # two groups of size two each
df <- c(1, 4) # dof for each of the two groups
U <- rgStudentcopula(n, groupings = groupings, df = df, scale = P)
(fit <- fitgStudentcopula(u = U, groupings = groupings, verbose = FALSE))
}
\keyword{distribution}