\name{riskmeasures}
\alias{VaR_nvmix}
\alias{ES_nvmix}
\title{Risk measures for normal variance mixtures}
\description{
  Estimation of value-at-risk and expected shortfall for univariate
  normal variance mixtures
}
\usage{
VaR_nvmix(level, qmix, loc = 0, scale = 1, control = list(), verbose = TRUE, ...)
ES_nvmix(level, qmix, loc = 0, scale = 1, control = list(), verbose = TRUE, ...)
}
\arguments{
  \item{level}{\eqn{n}-\code{\link{vector}} of confidence levels.}
  \item{qmix}{see \code{\link{pnvmix}()}.}
  \item{loc}{\code{\link{numeric}} location, see also \code{\link{pnvmix}()}}
  \item{scale}{\code{\link{numeric}} scale, see also \code{\link{pnvmix}()}}
  \item{control}{\code{\link{list}} specifying algorithm specific
    parameters; see \code{\link{get_set_param}()}.}
  \item{verbose}{
    \code{\link{logical}}
    indicating whether a warning is given if the required precision
    has not been reached.}
  \item{\dots}{additional arguments (for example, parameters) passed to
    the underlying mixing distribution when \code{qmix} is a
    \code{\link{character}} string or \code{\link{function}}, see also
    \code{\link{pnvmix}()}}
}
\value{
  \code{VaR_nvmix()} and \code{ES_nvmix()} return
  a \code{\link{numeric}} \eqn{n}-vector with the computed
  risk measures and in case of \code{ES_nvmix()} corresponding attributes
  \code{"abs. error"} and \code{"rel. error"}(error estimates of the RQMC estimator)
  and \code{"numiter"} (number of iterations).
}
\details{
  \code{VaR_nvmix} calls \code{\link{qnvmix}()}.

  The function \code{ES_nvmix()} estimates the expected shortfall using a
  randomized quasi Monte Carlo procedure by sampling from the mixing variable
  specified via \code{qmix} and and using the identity
  \eqn{\int_k^{\infty} x\phi(x)dx=\phi(k)} where \eqn{\phi(x)} denotes the
  density of a standard normal distribution.
  Algorithm specific paramaters (such as tolerances) can be conveniently passed
  via the \code{control} argument, see \code{\link{get_set_param}()} for more
  details.
}
\author{Erik Hintz, Marius Hofert and Christiane Lemieux}
\references{
  Hintz, E., Hofert, M. and Lemieux, C. (2021),
  Normal variance mixtures: Distribution, density and parameter estimation.
  \emph{Computational Statistics and Data Analysis} 157C, 107175.

  Hintz, E., Hofert, M. and Lemieux, C. (2022),
  Multivariate Normal Variance Mixtures in \R: The \R Package nvmix.
  \emph{Journal of Statistical Software}, \doi{10.18637/jss.v102.i02}.
}
\seealso{
  \code{\link{dnvmix}()}, \code{\link{pnvmix}()}, \code{\link{qnvmix}()},
  \code{\link{rnvmix}()}, \code{\link{get_set_param}()}
}
\examples{
## Example for inverse-gamma mixture (resulting in a t distribution) for
## which the expected shortfall admits a closed formula
set.seed(42) # reproducibility
level <- seq(from = 0.9, to = 0.95, by = 0.01)
df <- 4
## If 'qmix' is provided as string, ES_nvmix() uses the closed formula
ES1 <- ES_nvmix(level, qmix = "inverse.gamma", df = df)
## If 'qmix' is provided as function, the expected shortfall is estimated
ES2 <- ES_nvmix(level, qmix = function(u, df) 1/qgamma(1-u, shape = df/2, rate = df/2),
       df = df)
stopifnot(all.equal(ES1, ES2, tol = 1e-2, check.attributes = FALSE))
}
\keyword{distribution}