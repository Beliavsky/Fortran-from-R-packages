\name{skewstudent}
\alias{rskewt}
\alias{dskewt}
\alias{rskewtcopula}
\alias{dskewtcopula}
\title{Functionalities for the skew-t distribution and copula}
\description{
Sampling and density evaluation for the multivariate skew-t distribution and copula.
}
\usage{
rskewt(n, loc = rep(0, d), scale = diag(2), factor = NULL, gamma = rep(0, d),
       df = Inf, method = c("PRNG", "sobol", "ghalton"), skip = 0)
dskewt(x, loc = rep(0, d), scale = diag(2), gamma = rep(0, d), df,
       log = FALSE, scale.inv, ldet)

rskewtcopula(n, scale = diag(2), factor = NULL, gamma = rep(0, d), df = Inf,
             pseudo = TRUE, method = c("PRNG", "sobol", "ghalton"), skip = 0)
dskewtcopula(u, scale = diag(2), gamma = rep(0, d), df, log = FALSE,
             scale.inv, ldet)
}
\arguments{
  \item{u}{\eqn{(n, d)}-\code{\link{matrix}} of evaluation points or data;
    Have to be in \eqn{(0,1)}.}
  \item{x}{\eqn{(n, d)}-\code{\link{matrix}} of evaluation points or data}
  \item{n}{sample size \eqn{n} (positive integer).}
  \item{df}{positive degress of freedom; can also be \code{Inf} in which
    case the copula is interpreted as the Gaussian copula.}
  \item{loc}{location of length \eqn{d}.}
  \item{gamma}{Skewness-vector of dimension \eqn{d}; if \code{all(gamma == 0)},
  the classical t distribution or copula results.}
  \item{scale}{scale matrix (a covariance matrix entering the
    distribution as a parameter) of dimension \eqn{(d, d)} (defaults to
    \eqn{d = 2}). Note that \code{scale} must be positive definite,}
  \item{factor}{% more general here than for dnvmix() or pnvmix()
    \eqn{(d, d)}-\code{\link{matrix}} such that
    \code{factor \%*\% t(factor)} equals \code{scale}.
    If not provided, \code{factor} is internally determined via \code{\link{chol}()}}
  \item{scale.inv}{inverse of \code{scale}; if not provided, computed via
    \code{pd.solve(scale)}.}
  \item{ldet}{\code{log(det(scale))}; if not provided, computed via \code{pd.solve(scale)}.}
  \item{log}{\code{\link{logical}} indicating whether the logarithmic
    density is to be computed.}
  \item{pseudo}{\code{\link{logical}}; if \code{TRUE}, copula samples are computed
  via \code{pobs()} from a multivariate skew-t sample. If \code{FALSE}, the univariate skew t distribution functions are internally approximated via
  \code{integrate()}; see details below.}
  \item{method}{see \code{\link{rnvmix}()}.}
  \item{skip}{see \code{\link{rnvmix}()}.}
}
\value{
  \eqn{n}-vector of (log-)density values and \eqn{(n, d)}-matrix of samples, respectively.
}
\details{
Functionalities for sampling from the multivariate skew-\emph{t} distribution
and copula; the former has stochastic representation \eqn{\mu + W\gamma + \sqrt{W}AZ}
where \eqn{AA^T=scale}, \eqn{W} follows an inverse-gamma distrubution with
parameters \code{df/2} and is independent of the \eqn{d}-dimensional vector \eqn{Z}
following a standard multivariate normal distribution. When \eqn{gamma} is the
null-vector, the distribution becomes the multivariate \eqn{t} distribution.

A major computational challenge when working with the skew \emph{t} copula is
the lack of an available distribution and quantile function of the univariate
skew \emph{t} distribution. These are required in \code{rskewtcopula(, pobs = FALSE)}
and in \code{dskewtcopula()}. The unviarate skew \emph{t} distribution and
quantile functions are currently implemented as described Yoshiba, T. (2018).

The functions described here are currently being further developed to improve
stability, accuracy and speed, so that arguments may change in subsequent
versions of \code{nvmix}.
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

  Yoshiba, T. (2018).
  Maximum Likelihood Estimation of Skew-\emph{t} Copulas with Its Applications to Stock
  Returns.
  \emph{Journal of Statistical Computation and Simulation} 88 (13): 2489–2506.
}
\seealso{
  \code{\link{rStudent}()}, \code{\link{dStudent}()}, \code{\link{rStudentcopula}()}, \code{\link{dStudentcopula}()}}
\examples{
## Sampling from the skew-t copula

n <- 100 # sample size
d <- 10 # dimension
rho <- 0.5
scale <- matrix(rho, ncol = d, nrow = d)
diag(scale) <- 1 # scale
gamma <- rep(1, d) # skewness
df <- 7 # degrees-of-freedom parameter
set.seed(1) # same random numbers for both runs
system.time(samplecop_pobs <- rskewtcopula(n, scale = scale, gamma = gamma,
                                           df = df, pseudo = TRUE))
set.seed(1)
system.time(samplecop_pskewt <- rskewtcopula(n, scale = scale, gamma = gamma,
                                             df = df, pseudo = FALSE))
## Plot first two coordinates
layout(rbind(1:2))
plot(samplecop_pobs, xlab = expression(U[1]), ylab = expression(U[2]))
mtext("pobs = TRUE")
plot(samplecop_pskewt, xlab = expression(U[1]), ylab = expression(U[2]))
mtext("pobs = FALSE")
layout(1)
}
\keyword{distribution}