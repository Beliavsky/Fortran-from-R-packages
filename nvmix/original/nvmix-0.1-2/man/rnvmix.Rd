\name{rnvmix}
\alias{rnvmix}
\alias{rStudent}
\alias{rNorm}
\alias{rNorm_sumconstr}
\title{(Quasi-)Random Number Generation for Multivariate Normal Variance Mixtures}
\description{
  Generate vectors of random variates from multivariate normal variance
  mixtures (including Student \emph{t} and normal distributions).
}
\usage{
rnvmix(n, rmix, qmix, loc = rep(0, d), scale = diag(2),
       factor = NULL, method = c("PRNG", "sobol", "ghalton"),
       skip = 0, ...)

rStudent(n, df, loc = rep(0, d), scale = diag(2), factor = NULL,
         method = c("PRNG", "sobol", "ghalton"), skip = 0)
rNorm(n, loc = rep(0, d), scale = diag(2), factor = NULL,
      method = c("PRNG", "sobol", "ghalton"), skip = 0)
rNorm_sumconstr(n, weights, s, method = c("PRNG", "sobol", "ghalton"), skip = 0)
}
\arguments{
  \item{n}{sample size \eqn{n} (positive integer).}
  \item{rmix}{specification of the mixing variable \eqn{W}, see McNeil
    et al. (2015, Chapter 6) and Hintz et al. (2020),
    via a random number generator.
    This argument is ignored for \code{method = "sobol"} and
    \code{method = "ghalton"}. Supported are the
    following types of specification (see also the examples below):
    \describe{
      \item{\code{\link{character}}:}{\code{\link{character}} string
	specifying a supported distribution; currently available are
        \code{"constant"} (in which case \eqn{W = 1} and thus a sample
	from the multivariate normal distribution with mean vector
	\code{loc} and covariance matrix \code{scale} results) and
	\code{"inverse.gamma"} (in which case \eqn{W} is
	inverse gamma distributed with shape and rate parameters
	\code{df}/2 and thus the multivariate
	Student \emph{t} distribution with \code{df} degrees of freedom
	(required to be provided via the ellipsis argument) results).}
      \item{\code{\link{list}}:}{\code{\link{list}} of length at least
	one, where the first component is a \code{\link{character}}
	string specifying the base name of a distribution which can be
	sampled via prefix \code{"r"}; an example is \code{"exp"}
        for which \code{"rexp"} exists for sampling. If the list is
        of length larger than one, the remaining elements contain
        additional parameters of the distribution; for \code{"exp"},
	for example, this can be the parameter \code{rate}.}
      \item{\code{\link{function}}:}{\code{\link{function}}
	interpreted as a random number generator of the mixing
	variable \eqn{W}; additional arguments (such as parameters)
        can be passed via the ellipsis argument.}
      \item{\code{\link{numeric}}:}{\code{\link{numeric}} vector
	of length \code{n} providing a random sample of the mixing
	variable \eqn{W}.}
    }
  }
  \item{qmix}{specification of the mixing variable \eqn{W} via a quantile
    function. This argument is required
    for \code{method = "sobol"} and \code{method = "ghalton"}. Supported are
    the following types of specification (see also the examples below):
    \describe{
      \item{\code{\link{character}}:}{\code{\link{character}} string
	specifying a supported distribution; currently available are
        \code{"constant"} (in which case \eqn{W = 1} and thus a sample
	from the multivariate normal distribution with mean vector
	\code{loc} and covariance matrix \code{scale} results) and
	\code{"inverse.gamma"} (in which case \eqn{W} is
	inverse gamma distributed with shape and rate parameters
	\code{df}/2 and thus the multivariate
	Student \emph{t} distribution with \code{df} degrees of freedom
	(required to be provided via the ellipsis argument) results).}
      \item{\code{\link{list}}:}{\code{\link{list}} of length at least
	one, where the first component is a \code{\link{character}}
	string specifying the base name of a distribution which can be
	sampled via prefix \code{"q"}; an example is \code{"exp"}
        for which \code{"qexp"} exists for sampling. If the list is
        of length larger than one, the remaining elements contain
        additional parameters of the distribution; for \code{"exp"},
	for example, this can be the parameter \code{rate}.}
      \item{\code{\link{function}}:}{\code{\link{function}}
	interpreted as the quantile function of the mixing
	variable \eqn{W}; internally, sampling is then done with the
	inversion method by applying the provided function to U(0,1)
	random variates.}
    }
  }
  \item{df}{positive degress of freedom; can also be \code{Inf} in which
    case the distribution is interpreted as the multivariate normal
    distribution with mean vector \code{loc} and covariance matrix
    \code{scale}).}
  \item{loc}{location vector of dimension \eqn{d}; this equals the mean
    vector of a random vector following the specified normal variance
    mixture distribution if and only if the latter exists.}
  \item{scale}{scale matrix (a covariance matrix entering the
    distribution as a parameter) of dimension \eqn{(d, d)} (defaults to
    \eqn{d = 2});
    this equals the covariance matrix of a random vector following
    the specified normal variance mixture distribution divided by
    the expecation of the mixing variable \eqn{W} if and only if the
    former exists.
    Note that \code{scale} must be positive definite; sampling from
    singular normal variance mixtures can be achieved by providing
    \code{factor}.}
  \item{factor}{% more general here than for dnvmix() or pnvmix()
    \eqn{(d, k)}-\code{\link{matrix}} such that
    \code{factor \%*\% t(factor)} equals \code{scale}; the non-square
    case \eqn{k \neq d}{k != d} can be used to sample from singular
    normal variance mixtures. Note that this notation coincides with
    McNeil et al. (2015, Chapter 6). If not provided, \code{factor} is
    internally determined via \code{\link{chol}()} (and multiplied
    from the right to an \eqn{(n, k)}-matrix of independent
    standard normals to obtain a sample from a multivariate normal
    with zero mean vector and covariance matrix \code{scale}).}
  \item{method}{\code{\link{character}} string indicating the method
    to be used to obtain the sample. Available are:
    \describe{
      \item{\code{"PRNG"}:}{pseudo-random numbers,}
      \item{\code{"sobol"}:}{Sobol' sequence,}
      \item{\code{"ghalton"}:}{generalized Halton sequence.}
    }
    If \code{method = "PRNG"}, either \code{qmix} or \code{rmix} can be
    provided. If both are provided, \code{rmix} is used and \code{qmix}
    ignored. For the other two methods, sampling is done via inversion,
    hence \code{qmix} has to be provided and \code{rmix} is ignored.
  }
  \item{skip}{\code{\link{integer}} specifying the number of points
    to be skipped when \code{method = "sobol"}, see also example below.
  }
  \item{weights}{\eqn{d}-\code{\link{numeric}} vector of weights.}
  \item{s}{\code{\link{numeric}} vector of length 1 or \code{n} giving the
  value of the constrained sum; see below under details.}
  \item{\dots}{additional arguments (for example, parameters) passed to
    the underlying mixing distribution when \code{rmix} or \code{qmix}
    is a \code{\link{character}} string or \code{\link{function}}.}
}
\value{
  \code{rnvmix()} returns an \eqn{(n, d)}-\code{\link{matrix}}
  containing \eqn{n} samples of the specified (via \code{mix})
  \eqn{d}-dimensional multivariate normal variance mixture with
  location vector \code{loc} and scale matrix \code{scale}
  (a covariance matrix).

  \code{rStudent()} returns samples from the \eqn{d}-dimensional
  multivariate Student \emph{t} distribution with location vector
  \code{loc} and scale matrix \code{scale}.

  \code{rNorm()} returns samples from the \eqn{d}-dimensional
  multivariate normal distribution with mean vector
  \code{loc} and covariance matrix \code{scale}.

  \code{rNorm_sumconstr()} returns samples from the \eqn{d}-dimensional
  multivariate normal distribution conditional on the weighted sum being
  constrained to \code{s}.
}
\details{
  Internally used is \code{factor}, so \code{scale} is not required
  to be provided if \code{factor} is given.

  The default factorization used to obtain \code{factor} is the Cholesky
  decomposition via \code{\link{chol}()}. To this end, \code{scale}
  needs to have full rank.

  Sampling from a singular normal variance mixture distribution can be
  achieved by providing \code{factor}.

  The number of rows of \code{factor} equals the dimension \eqn{d} of
  the sample. Typically (but not necessarily), \code{factor} is square.

  \code{rStudent()} and \code{rNorm()} are wrappers of
  \code{rnvmix(, qmix = "inverse.gamma", df = df)} and
  \code{rnvmix(, qmix = "constant", df = df)}, respectively.

  The function \code{rNorm_sumconstr()} can be used to sample from the
  multivariate standard normal distribution under a weighted sum constraint;
  the implementation is based on Algorithm 1 in Vrins (2018). Let
  \eqn{Z = (Z_1,\dots,Z_d)~N_d(0, I_d)}. The function \code{rNorm_sumconstr()}
  then samples from \eqn{Z | w^T Z = s} where \eqn{w} and \eqn{s} correspond
  to the arguments \code{weights} and \code{s}. If supplied \code{s} is a vector
  of length \code{n}, the i'th row of the returned matrix uses the constraint
  \eqn{w^T Z = s_i} where \eqn{s_i} is the i'th element in \code{s}.

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

  Vrins, E. (2018)
  Sampling the Multivariate Standard Normal Distribution under a Weighted
  Sum Constraint.
  \emph{Risks} 6(3), 64.
}
\seealso{
  \code{\link{dnvmix}()}, \code{\link{pnvmix}()}
}
\examples{
### Examples for rnvmix() ######################################################

## Generate a random correlation matrix in d dimensions
d <- 3
set.seed(157)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A \%*\% t(A))

## Draw random variates and compare
df <- 3.5
n <- 1000
set.seed(271)
X  <- rnvmix(n, rmix = "inverse.gamma", df = df, scale = P) # with scale
set.seed(271)
X. <- rnvmix(n, rmix = "inverse.gamma", df = df, factor = t(chol(P))) # with factor
stopifnot(all.equal(X, X.))

## Checking df = Inf
set.seed(271)
X <- rnvmix(n, rmix = "constant", scale = P) # normal
set.seed(271)
X. <- rnvmix(n, rmix = "inverse.gamma", scale = P, df = Inf) # t_infinity
stopifnot(all.equal(X, X.))

## Univariate case (dimension = number of rows of 'factor' = 1 here)
set.seed(271)
X.1d  <- rnvmix(n, rmix = "inverse.gamma", df = df, factor = 1/2)
set.seed(271)
X.1d. <- rnvmix(n, rmix = "inverse.gamma", df = df, factor = 1)/2 # manual scaling
stopifnot(all.equal(X.1d, X.1d.))

## Checking different ways of providing 'mix'
## 1) By providing a character string (and corresponding ellipsis arguments)
set.seed(271)
X.mix1 <- rnvmix(n, rmix = "inverse.gamma", df = df, scale = P)
## 2) By providing a list; the first element has to be an existing distribution
##    with random number generator available with prefix "r"
rinverse.gamma <- function(n, df) 1 / rgamma(n, shape = df/2, rate = df/2)
set.seed(271)
X.mix2 <- rnvmix(n, rmix = list("inverse.gamma", df = df), scale = P)
## 3) The same without extra arguments (need the extra list() here to
##    distinguish from Case 1))
rinverseGamma <- function(n) 1 / rgamma(n, shape = df/2, rate = df/2)
set.seed(271)
X.mix3 <- rnvmix(n, rmix = list("inverseGamma"), scale = P)
## 4) By providing a quantile function
##    Note: P(1/Y <= x) = P(Y >= 1/x) = 1-F_Y(1/x) = y <=> x = 1/F_Y^-(1-y)
set.seed(271)
X.mix4 <- rnvmix(n, qmix = function(p) 1/qgamma(1-p, shape = df/2, rate = df/2),
                 scale = P)
## 5) By providing random variates
set.seed(271) # if seed is set here, results are comparable to the above methods
W <- rinverse.gamma(n, df = df)
X.mix5 <- rnvmix(n, rmix = W, scale = P)
## Compare (note that X.mix4 is not 'all equal' with X.mix1 or the other samples)
## since rgamma() != qgamma(runif()) (or qgamma(1-runif()))
stopifnot(all.equal(X.mix2, X.mix1),
          all.equal(X.mix3, X.mix1),
          all.equal(X.mix5, X.mix1))

## For a singular normal variance mixture:
## Need to provide 'factor'
A <- matrix( c(1, 0, 0, 1, 0, 1), ncol = 2, byrow = TRUE)
stopifnot(all.equal(dim(rnvmix(n, rmix = "constant", factor = A)),    c(n, 3)))
stopifnot(all.equal(dim(rnvmix(n, rmix = "constant", factor = t(A))), c(n, 2)))

## Using 'skip'. Need to reset the seed everytime to get the same shifts in "sobol".
## Note that when using method = "sobol", we have to provide 'qmix' instead of 'rmix'.
set.seed(271)
X.skip0 <- rnvmix(n, qmix = "inverse.gamma", df = df, scale = P, method = "sobol")
set.seed(271)
X.skip1 <- rnvmix(n, qmix = "inverse.gamma", df = df, scale = P, method = "sobol",
                  skip = n)
set.seed(271)
X.wo.skip <- rnvmix(2*n, qmix = "inverse.gamma", df = df, scale = P, method = "sobol")
X.skip <- rbind(X.skip0, X.skip1)
stopifnot(all.equal(X.wo.skip, X.skip))


### Examples for rStudent() and rNorm() ########################################

## Draw N(0, P) random variates by providing scale or factor and compare
n <- 1000
set.seed(271)
X.n  <- rNorm(n, scale = P) # providing scale
set.seed(271)
X.n. <- rNorm(n, factor = t(chol(P))) # providing the factor
stopifnot(all.equal(X.n, X.n.))

## Univariate case (dimension = number of rows of 'factor' = 1 here)
set.seed(271)
X.n.1d  <- rNorm(n, factor = 1/2)
set.seed(271)
X.n.1d. <- rNorm(n, factor = 1)/2 # manual scaling
stopifnot(all.equal(X.n.1d, X.n.1d.))

## Draw t_3.5 random variates by providing scale or factor and compare
df <- 3.5
n <- 1000
set.seed(271)
X.t  <- rStudent(n, df = df, scale = P) # providing scale
set.seed(271)
X.t. <- rStudent(n, df = df, factor = t(chol(P))) # providing the factor
stopifnot(all.equal(X.t, X.t.))

## Univariate case (dimension = number of rows of 'factor' = 1 here)
set.seed(271)
X.t.1d  <- rStudent(n, df = df, factor = 1/2)
set.seed(271)
X.t.1d. <- rStudent(n, df = df, factor = 1)/2 # manual scaling
stopifnot(all.equal(X.t.1d, X.t.1d.))

## Check df = Inf
set.seed(271)
X.t <- rStudent(n, df = Inf, scale = P)
set.seed(271)
X.n <- rNorm(n, scale = P)
stopifnot(all.equal(X.t, X.n))

### Examples for rNorm_sumconstr() #############################################
set.seed(271)
weights <- c(1, 1)
Z.constr <- rNorm_sumconstr(n, weights = c(1, 1), s = 2)
stopifnot(all(rowSums(Z.constr ) == 2))
plot(Z.constr , xlab = expression(Z[1]), ylab = expression(Z[2]))
}
\keyword{distribution}