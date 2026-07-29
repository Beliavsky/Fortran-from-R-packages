\name{dgnvmix}
\alias{dgnvmix}
\alias{dgStudent}
\title{Density of Grouped Normal Variance Mixtures}
\description{
  Evaluating grouped normal variance mixture density functions
  (including Student \emph{t} with multiple degrees-of-freedom).
}
\usage{
dgnvmix(x, groupings = 1:d, qmix, loc = rep(0, d), scale = diag(d), factor = NULL,
        factor.inv = NULL, control = list(), log = FALSE, verbose = TRUE, ...)
dgStudent(x, groupings = 1:d, df, loc = rep(0, d), scale = diag(d), factor = NULL,
          factor.inv = NULL, control = list(), log = FALSE, verbose = TRUE)
}
\arguments{
  \item{x}{see \code{\link{dnvmix}()}.}
  \item{groupings}{see \code{\link{pgnvmix}()}.}
  \item{qmix}{specification of the mixing variables \eqn{W_i} via quantile
  functions; see \code{\link{pgnvmix}()}.}
  \item{loc}{see \code{\link{pnvmix}()}.}
  \item{scale}{see \code{\link{pnvmix}()}; must be positive definite.}
  \item{factor}{\eqn{(d, d)} lower triangular \code{\link{matrix}}
    such that \code{factor \%*\% t(factor)} equals \code{scale}. Internally
    used is \code{factor.inv}.}
  \item{factor.inv}{inverse of \code{factor}; if not provided, computed via
    \code{solve(factor)}.}
  \item{df}{\code{\link{vector}} of length \code{length(unique(groupings))} so that
  variable \code{i} has degrees-of-freedom \code{df[groupings[i]]}; all elements
  must be positive
  and can be \code{Inf}, in which case the corresponding marginal is normally
  distributed.}
  \item{control}{\code{\link{list}} specifying algorithm specific
    parameters; see \code{\link{get_set_param}()}.}
  \item{log}{\code{\link{logical}} indicating whether the logarithmic
    density is to be computed.}
  \item{verbose}{see \code{\link{pnvmix}()}.}
  \item{\dots}{additional arguments (for example, parameters) passed
    to the underlying mixing distribution when \code{qmix} is a
    \code{\link{character}} string or an element of \code{qmix} is a
    \code{\link{function}}.}
}
\value{
  \code{dgnvmix()} and \code{dgStudent()} return a
  \code{\link{numeric}} \eqn{n}-vector with the computed density values
  and corresponding attributes \code{"abs. error"} and \code{"rel. error"}
  (error estimates of the RQMC estimator) and \code{"numiter"} (number of iterations).
}
\details{
  Internally used is \code{factor.inv}, so \code{factor} and \code{scale} are not
  required to be provided (but allowed for consistency with other functions in the
  package).

  \code{dgStudent()} is a wrapper of
  \code{dgnvmix(, qmix = "inverse.gamma", df = df)}. If there is no grouping,
  the analytical formula for the density of a multivariate \emph{t} distribution
  is used.

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
}
\seealso{
  \code{\link{rgnvmix}()}, \code{\link{pgnvmix}()}, \code{\link{get_set_param}()}
}
\examples{
n <- 100 # sample size to generate evaluation points

### 1. Inverse-gamma mixture
## 1.1. Grouped t with mutliple dof

d <- 3 # dimension
set.seed(157)
A <- matrix(runif(d * d), ncol = d)
P <- cov2cor(A \%*\% t(A)) # random scale matrix
df <- c(1.1, 2.4, 4.9) # dof for margin i
groupings <- 1:d
x <- rgStudent(n, df = df, scale = P) # evaluation points for the density

### Call 'dgnvmix' with 'qmix' a string:
set.seed(12)
dgt1 <- dgnvmix(x, qmix = "inverse.gamma", df = df, scale = P)
### Version providing quantile functions of the mixing distributions as list
qmix_ <- function(u, df) 1 / qgamma(1-u, shape = df/2, rate = df/2)
qmix <- list(function(u) qmix_(u, df = df[1]), function(u) qmix_(u, df = df[2]),
             function(u) qmix_(u, df = df[3]))
set.seed(12)
dgt2 <- dgnvmix(x, groupings = groupings, qmix = qmix, scale = P)
### Similar, but using ellipsis argument:
qmix <- list(function(u, df1) qmix_(u, df1), function(u, df2) qmix_(u, df2),
             function(u, df3) qmix_(u, df3))
set.seed(12)
dgt3 <- dgnvmix(x, groupings = groupings, qmix = qmix, scale = P, df1 = df[1],
                df2 = df[2], df3 = df[3])
### Using the wrapper 'dgStudent()'
set.seed(12)
dgt4 <- dgStudent(x, groupings = groupings, df = df, scale = P)
stopifnot(all.equal(dgt1, dgt2, tol = 1e-5, check.attributes = FALSE),
          all.equal(dgt1, dgt3, tol = 1e-5, check.attributes = FALSE),
          all.equal(dgt1, dgt4, tol = 1e-5, check.attributes = FALSE))


## 1.2 Classical multivariate t

df <- 2.4
groupings <- rep(1, d) # same df for all components
x <- rStudent(n, scale = P, df = df) # evaluation points for the density
dt1 <- dStudent(x, scale = P, df = df, log = TRUE) # uses analytical formula
## If 'qmix' provided as string and no grouping, dgnvmix() uses analytical formula
dt2 <- dgnvmix(x, qmix = "inverse.gamma", groupings = groupings, df = df, scale = P, log = TRUE)
stopifnot(all.equal(dt1, dt2))
## Provide 'qmix' as a function to force estimation in 'dgnvmix()'
dt3 <- dgnvmix(x, qmix = qmix_, groupings = groupings, df = df, scale = P, log = TRUE)
stopifnot(all.equal(dt1, dt3, tol = 5e-4, check.attributes = FALSE))

### 2. More complicated mixutre

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
x <- rgnvmix(n, groupings = groupings, qmix = qmix, scale = P)
dg <- dgnvmix(x, groupings = groupings, qmix = qmix, scale = P, log = TRUE)
}
\keyword{distribution}