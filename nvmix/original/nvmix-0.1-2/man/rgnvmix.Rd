\name{rgnvmix}
\alias{rgnvmix}
\alias{rgStudent}
\title{(Quasi-)Random Number Generator for Grouped Normal Variance Mixtures}
\description{
  Generate vectors of random variates from grouped normal variance
  mixtures (including Student \emph{t} with multiple degrees-of-freedom).
}
\usage{
rgnvmix(n, qmix, groupings = 1:d, loc = rep(0, d), scale = diag(2),
       factor = NULL, method = c("PRNG", "sobol", "ghalton"), skip = 0, ...)
rgStudent(n, groupings = 1:d, df, loc = rep(0, d), scale = diag(2),
          factor = NULL, method = c("PRNG", "sobol", "ghalton"), skip = 0)
}
\arguments{
  \item{n}{sample size \eqn{n} (positive integer).}
  \item{qmix}{specification of the mixing variables \eqn{W_i}; see
    \code{\link{pgnvmix}()}.}
  \item{groupings}{\code{\link{vector}} specifying the group structure; see
  \code{\link{pgnvmix}()}.}
  \item{df}{\code{\link{vector}} specifying the degrees-of-freedom; see
   see \code{\link{pgStudent}()}.}
  \item{loc}{see \code{\link{pgnvmix}()}.}
  \item{scale}{see \code{\link{pgnvmix}()}. \code{scale} must be positive definite;
  sampling from singular normal variance mixtures can be achieved by providing
    \code{factor}.}
  \item{factor}{% more general here than for dnvmix() or pnvmix()
    see \code{\link{rnvmix}()}.}
  \item{method}{see \code{\link{rnvmix}()}.}
  \item{skip}{see \code{\link{rnvmix}()}.}
  \item{\dots}{additional arguments (for example, parameters) passed to
    the underlying mixing distribution when \code{qmix}
    is a \code{\link{character}} string or an element of \code{qmix} is
    a \code{\link{function}}.}
}
\value{
  \code{rgnvmix()} returns an \eqn{(n, d)}-\code{\link{matrix}}
  containing \eqn{n} samples of the specified (via \code{qmix})
  \eqn{d}-dimensional grouped normal variance mixture with
  location vector \code{loc} and scale matrix \code{scale}
  (a covariance matrix).

  \code{rgStudent()} returns samples from the \eqn{d}-dimensional
  multivariate \emph{t} distribution with multiple degrees-of-freedom
  specified by \code{df}, location vector
  \code{loc} and scale matrix \code{scale}.
}
\details{
  Internally used is \code{factor}, so \code{scale} is not required
  to be provided if \code{factor} is given.

  The default factorization used to obtain \code{factor} is the Cholesky
  decomposition via \code{\link{chol}()}. To this end, \code{scale}
  needs to have full rank.

  \code{rgStudent()} is a wrapper of
  \code{rgnvmix(, qmix = "inverse.gamma", df = df)}.
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
  \code{\link{rnvmix}()}, \code{\link{pgnvmix}()}
}
\examples{
n <- 1000 # sample size

## Generate a random correlation matrix in d dimensions
d <- 2
set.seed(157)
A <- matrix(runif(d * d), ncol = d)
scale <- cov2cor(A \%*\% t(A))

## Example 1: Exponential mixture
## Let W_1 ~ Exp(1), W_2 ~ Exp(10)
rates <- c(1, 10)
#qmix <- list(list("exp", rate = rates[1]), list("exp", rate = rates[2]))
qmix <- lapply(1:2, function(i) list("exp", rate = rates[i]))
set.seed(1)
X.exp1 <- rgnvmix(n, qmix = qmix, scale = scale)
## For comparison, consider NVM distribution with W ~ Exp(1)
set.seed(1)
X.exp2 <- rnvmix(n, qmix = list("exp", rate = rates[1]), scale = scale)
## Plot both samples with the same axes
opar <- par(no.readonly = TRUE)
par(mfrow=c(1,2))
plot(X.exp1, xlim = range(X.exp1, X.exp2), ylim = range(X.exp1, X.exp2),
     xlab = expression(X[1]), ylab = expression(X[2]))
mtext("Two groups with rates 1 and 10")
plot(X.exp2, xlim = range(X.exp1, X.exp2), ylim = range(X.exp1, X.exp2),
     xlab = expression(X[1]), ylab = expression(X[2]))
mtext("One group with rate 1")
par(opar)

## Example 2: Exponential + Inverse-gamma mixture
## Let W_1 ~ Exp(1), W_2 ~ IG(1.5, 1.5) (=> X_2 ~ t_3 marginally)
df <- 3
qmix <- list(list("exp", rate = rates[1]),
             function(u, df) 1/qgamma(1-u, shape = df/2, rate = df/2))
set.seed(1)
X.mix1 <- rgnvmix(n, qmix = qmix, scale = scale, df = df)
plot(X.mix1, xlab = expression(X[1]), ylab = expression(X[2]))

## Example 3: Mixtures in d > 2
d <- 5
set.seed(157)
A <- matrix(runif(d * d), ncol = d)
scale <- cov2cor(A \%*\% t(A))

## Example 3.1: W_i ~ Exp(i), i = 1,...,d
qmix <- lapply(1:d, function(i) list("exp", rate = i))
set.seed(1)
X.mix2 <- rgnvmix(n, qmix = qmix, scale = scale)

## Example 3.2: W_1, W_2 ~ Exp(1), W_3, W_4, W_5 ~ Exp(2)
## => 2 groups, so we need two elements in 'qmix'
qmix <- lapply(1:2, function(i) list("exp", rate = i))
groupings <- c(1, 1, 2, 2, 2)
set.seed(1)
X.mix3 <- rgnvmix(n, qmix = qmix, groupings = groupings, scale = scale)

## Example 3.3: W_1, W_3 ~ IG(1, 1), W_2, W_4 ~ IG(2, 2), W_5 = 1
## => X_1, X_3 ~ t_2; X_2, X_4 ~ t_4, X_5 ~ N(0, 1)
qmix <- list(function(u, df1) 1/qgamma(1-u, shape = df1/2, rate = df1/2),
             function(u, df2) 1/qgamma(1-u, shape = df2/2, rate = df2/2),
             function(u) rep(1, length(u)))
groupings = c(1, 2, 1, 2, 3)
df = c(2, 4, Inf)
set.seed(1)
X.t1 <- rgnvmix(n, qmix = qmix, groupings = groupings, scale = scale,
                df1 = df[1], df2 = df[2])

## This is equivalent to calling 'rgnmvix' with 'qmix = "inverse.gamma"'
set.seed(1)
X.t2 <- rgnvmix(n, qmix = "inverse.gamma", groupings = groupings, scale = scale,
                df = df)

## Alternatively, one can use the user friendly wrapper 'rgStudent()'
set.seed(1)
X.t3 <- rgStudent(n, df = df, groupings = groupings, scale = scale)

stopifnot(all.equal(X.t1, X.t2), all.equal(X.t1, X.t3))
}
\keyword{distribution}