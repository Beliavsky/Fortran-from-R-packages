\name{pgnvmix}
\alias{pgnvmix}
\alias{pgStudent}
\title{Distribution Function of Grouped Multivariate Normal Variance Mixtures}
\description{
  Evaluating grouped and generalized multivariate normal variance mixture distribution
  functions (including Student \emph{t} with multiple degrees-of-freedom).
}
\usage{
pgnvmix(upper, lower = matrix(-Inf, nrow = n, ncol = d), groupings = 1:d, qmix,
        rmix, loc = rep(0, d), scale = diag(d), standardized = FALSE,
        control = list(), verbose = TRUE, ...)

pgStudent(upper, lower = matrix(-Inf, nrow = n, ncol = d), groupings = 1:d, df,
          loc = rep(0, d), scale = diag(d), standardized = FALSE,
          control = list(), verbose = TRUE)
}
\arguments{
  \item{upper}{see \code{\link{pnvmix}()}.}
  \item{lower}{see \code{\link{pnvmix}()}.}
  \item{groupings}{\eqn{d}-\code{\link{vector}}. Specification of the groupings
  so that variable \eqn{i} has mixing variable \eqn{W_k} where
  \code{k = groupings[i]}. If \code{groupings = 1:d}, each variable has a different
  mixing distribution. }
  \item{qmix}{specification of the mixing variables \eqn{W_i} via quantile
  functions; see McNeil et al. (2015, Chapter 6) and Hintz et al. (2020).
  Supported are the following types of
    specification (see also the examples below):
    \describe{\item{\code{\link{character}}:}{\code{\link{character}} string
	specifying a supported distribution; currently available are
	\code{"inverse.gamma"} (in which case \eqn{W_i} is
	inverse gamma distributed with shape and rate parameters
	\code{df[groupings[i]]}/2 and a multivariate
	Student \emph{t} distribution multiple degreess-of-freedom
	results) and
	\code{"pareto"} (in which case \eqn{W_i} is
	Pareto distributed with scale equal to unity and shape equal to
	\code{alpha[groupings[i]]}. \code{alpha} and \code{df} must be of length
	\code{length(unique(groupings))} and need to be provided via the ellipsis
        argument).}
      \item{\code{\link{list}}:}{\code{\link{list}} of length
      \code{length(unique(groupings))} (number of different mixing distributions).
      Element \eqn{i} of this list specifies the mixing variable for component
      \code{groupings[i]}. Each element of this list can be
      \describe{\item{\code{\link{list}}:}{a \code{\link{list}} of length at
      least one, where the first component is a \code{\link{character}}
      string specifying the base name of a distribution whose
	   quantile function can be accessed via the prefix \code{"q"}. An example
	   \code{"exp"} for which \code{"qexp"} exists.
	   If the list is of length larger than one, the remaining elements
	   contain additional parameters of the distribution; for \code{"exp"},
      for example, this can be the parameter \code{rate}.}
      \item{\code{\link{function}}:}{\code{\link{function}}
	   interpreted as the quantile function or random number generator of the mixing
	   variable \eqn{W_i}}}}
      }
   }
  \item{rmix}{only allowed when \code{groupings = rep(1, d)} in which case
  \code{pgnvmix()} is equivalent to \code{\link{pnvmix}()}; see \code{\link{pnvmix}()}.}
  \item{df}{\code{\link{vector}} of length \code{length(unique(groupings))} so that
  variable \code{i} has degrees-of-freedom \code{df[groupings[i]]}; all elements
  must be positive
  and can be \code{Inf}, in which case the corresponding marginal is normally
  distributed.}
  \item{loc}{see \code{\link{pnvmix}()}.}
  \item{scale}{see \code{\link{pnvmix}()}; must be positive definite.}
  \item{standardized}{see \code{\link{pnvmix}()}.}
  \item{control}{\code{\link{list}} specifying algorithm specific
    parameters; see \code{\link{get_set_param}()}.}
  \item{verbose}{see \code{\link{pnvmix}()}.}
  \item{\dots}{additional arguments (for example, parameters) passed
    to the underlying mixing distribution when \code{qmix} is a
    \code{\link{character}} string or an element of \code{qmix} is a
    \code{\link{function}}.}
}
\value{
  \code{pgnvmix()} and \code{pgStudent()} return a
  \code{\link{numeric}} \eqn{n}-vector with the computed probabilities
  and corresponding attributes \code{"abs. error"} and \code{"rel. error"}
  (error estimates of the RQMC estimator) and \code{"numiter"} (number of iterations).
}
\details{
  One should highlight that evaluating grouped normal variance mixtures is a
  non-trivial tasks which, at the time of development of \pkg{nvmix},
  was not available in \R before, not even the special case of a
  multivariate Student \emph{t} distribution for non-integer degrees of
  freedoms, which frequently appears in applications in finance,
  insurance and risk management after estimating such distributions.

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

  \code{pgStudent()} is a wrapper of
  \code{pgnvmix(, qmix = "inverse.gamma", df = df)}.
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
}
\seealso{
  \code{\link{rgnvmix}()}, \code{\link{dgnvmix}()}, \code{\link{get_set_param}()}
}
\examples{
### Examples for pgnvmix() #####################################################

## 1. Inverse-gamma mixture (=> distribution is grouped t with mutliple dof)
d <- 3
set.seed(157)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A \%*\% t(A))
a <- -3 * runif(d) * sqrt(d) # random lower limit
b <-  3 * runif(d) * sqrt(d) # random upper limit
df <- c(1.1, 2.4, 4.9) # dof for margin i
groupings <- 1:d

### Call 'pgnvmix' with 'qmix' a string:
set.seed(12)
(pgt1 <- pgnvmix(b, lower = a, groupings = groupings, qmix = "inverse.gamma",
                df = df, scale = P))
### Version providing quantile functions of the mixing distributions as list
qmix_ <- function(u, df) 1 / qgamma(1-u, shape = df/2, rate = df/2)
qmix <- list(function(u) qmix_(u, df = df[1]), function(u) qmix_(u, df = df[2]),
             function(u) qmix_(u, df = df[3]))
set.seed(12)
(pgt2 <- pgnvmix(b, lower = a, groupings = groupings, qmix = qmix, scale = P))
### Similar, but using ellipsis argument:
qmix <- list(function(u, df1) qmix_(u, df1), function(u, df2) qmix_(u, df2),
             function(u, df3) qmix_(u, df3))
set.seed(12)
(pgt3 <- pgnvmix(b, lower = a, groupings = groupings, qmix = qmix,
                        scale = P, df1 = df[1], df2 = df[2], df3 = df[3]))
## Version using the user friendly wrapper 'pgStudent()'
set.seed(12)
(pgt4 <- pgStudent(b, lower = a, groupings = groupings, scale = P, df = df))
stopifnot(all.equal(pgt1, pgt2, tol = 1e-4, check.attributes = FALSE),
          all.equal(pgt2, pgt3), all.equal(pgt1, pgt4))

## 2. More complicated mixutre
## Let W1 ~ IG(1, 1), W2 = 1, W3 ~ Exp(1), W4 ~ Par(2, 1), W5 = W1, all comonotone
## => X1 ~ t_2; X2 ~ normal; X3 ~ Exp-mixture; X4 ~ Par-mixture; X5 ~ t_2

d <- 5
set.seed(157)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A \%*\% t(A))
b <- 3 * runif(d) * sqrt(d) # random upper limit
groupings <- c(1, 2, 3, 4, 1) # since W_5 = W_1
qmix <- list(function(u) qmix_(u, df = 2), function(u) rep(1, length(u)),
             list("exp", rate=1), function(u) (1-u)^(-1/2)) # length 4 (# of groups)
pg1 <- pgnvmix(b, groupings = groupings, qmix = qmix, scale = P)
stopifnot(all.equal(pg1, 0.78711, tol = 5e-6, check.attributes = FALSE))
}
\keyword{distribution}