\name{qnvmix}
\alias{qnvmix}
\title{Quantile Function of a univariate Normal Variance Mixture Distribution}
\description{
  Evaluating multivariate normal variance mixture distribution functions
  (including normal and Student \emph{t} for non-integer degrees of freedom).
}
\usage{
qnvmix(u, qmix, control = list(),
       verbose = TRUE, q.only = TRUE, stored.values = NULL, ...)
}
\arguments{
  \item{u}{vector of probabilities .}
  \item{qmix}{specification of the mixing variable \eqn{W}; see
    \code{\link{pnvmix}()} for details and examples.}
  \item{control}{\code{\link{list}} specifying algorithm specific
    parameters; see \code{\link{get_set_param}()}.}
  \item{verbose}{\code{\link{logical}}, if \code{TRUE} a warning is printed if
    one of the error tolerances is not met.}
  \item{q.only}{\code{\link{logical}}. If \code{TRUE}, only the quantiles are
    returned; if \code{FALSE}, see Section 'value' below.}
  \item{stored.values}{\code{\link{matrix}} with 3 columns of the form
    \eqn{[x, F(x), logf(x)]} where \eqn{F()} and \eqn{logf()} are the distribution-
    and log-density function of the distribution specified in \code{qmix}.
    If provided it is used to determine starting values for internal newton proceudures.
    Only very basic checking is done.}
  \item{\dots}{additional arguments containing parameters of
    mixing distributions when \code{qmix} is a \code{\link{character}}
    string.}
}
\value{
  If \code{q.only = TRUE} a vector of the same length as \code{u} with
  entries \eqn{q_i} where \eqn{q_i} satisfies \eqn{q_i = inf_x { F(x)
  \ge u_i}} where \eqn{F(x)} the univariate df of the normal variance
  mixture specified via \code{qmix};

  if \code{q.only = FALSE} a list of four:
  \describe{
    \item{\code{$q}:}{Vector of quantiles,}
    \item{\code{$log.density}:}{vector log-density values at \code{q},}
    \item{\code{$computed.values}:}{matrix with 3 columns [x, F(x), logf(x)];
      see details above,}
    \item{\code{$newton.iterations}:}{vector giving the number of Newton
      iterations needed for \code{u[i]}.}
      }
}
\details{
  This function uses a Newton procedure to estimate the quantile of the
  specified univariate normal variance mixture distribution. Internally,
  a randomized quasi-Monte Carlo (RQMC) approach is used to estimate the
  distribution and (log)density function; the method is similar to the
  one in \code{\link{pnvmix}()} and \code{\link{dnvmix}()}. The result depends
  slightly on \code{.random.seed}.

  Internally, symmetry is used for \eqn{u \le 0.5}. Function values
  (i.e., df and log-density values) are stored and reused to get good
  starting values. These values are returned if \code{q.only = FALSE}
  and can be re-used by passing it to \code{qnvmix()} via the argument
  \code{stored.values}; this can significantly reduce run-time.

  Accuracy and run-time depend on both the magnitude of \eqn{u} and on
  how heavy the tail of the underlying distributions is.  Numerical
  instabilities can occur for values of \eqn{u} close to 0 or 1,
  especially when the tail of the distribution is heavy.

  If \code{q.only = FALSE} the log-density values of the underlying
  distribution evaluated at the estimated quantiles are returned as
  well: This can be useful for copula density evaluations where both
  quantities are needed.

  Underlying algorithm specific parameters can be changed via the \code{control}
  argument, see \code{\link{get_set_param}()} for details.
}
\seealso{
  \code{\link{dnvmix}()}, \code{\link{rnvmix}()}, \code{\link{pnvmix}()}
}
\author{Erik Hintz, Marius Hofert and Christiane Lemieux}
\references{
  Hintz, E., Hofert, M. and Lemieux, C. (2021),
  Normal variance mixtures: Distribution, density and parameter estimation.
  \emph{Computational Statistics and Data Analysis} 157C, 107175.

  Hintz, E., Hofert, M. and Lemieux, C. (2022),
  Multivariate Normal Variance Mixtures in \R: The \R Package nvmix.
  \emph{Journal of Statistical Software}, \doi{10.18637/jss.v102.i02}.

  McNeil, A. J., Frey, R., and Embrechts, P. (2015).
  \emph{Quantitative Risk Management: Concepts, Techniques, Tools}.
  Princeton University Press.
}
\examples{
## Evaluation points
u <- seq(from = 0.05, to = 0.95, by = 0.025)
set.seed(271) # for reproducibility

## Evaluate the t_{1.4} quantile function
df <- 1.4
qmix. <- function(u) 1/qgamma(1-u, shape = df/2, rate = df/2)
## If qmix = "inverse.gamma", qt() is being called
qt1 <- qnvmix(u, qmix = "inverse.gamma", df = df)
## Estimate quantiles (without using qt())
qt1. <- qnvmix(u, qmix = qmix., q.only = FALSE)
stopifnot(all.equal(qt1, qt1.$q, tolerance = 2.5e-3))
## Look at absolute error:
abs.error <- abs(qt1 - qt1.$q)
plot(u, abs.error, type = "l", xlab = "u", ylab = "Absolute error")
## Now do this again but provide qt1.$stored.values, in which case at most
## one Newton iteration will be needed:
qt2 <- qnvmix(u, qmix = qmix., stored.values = qt1.$computed.values, q.only = FALSE)
stopifnot(max(qt2$newton.iterations) <= 1)

}
\keyword{distribution}