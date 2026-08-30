\name{pmvnorm}
\alias{pmvnorm}
\title{ Multivariate Normal Distribution }
\description{

Computes the distribution function of the multivariate normal
distribution for arbitrary limits and correlation matrices.

}
\usage{
pmvnorm(lower=-Inf, upper=Inf, mean=rep(0, length(lower)),
           corr=NULL, sigma=NULL, algorithm = GenzBretz(), keepAttr=TRUE, 
           seed = NULL, \dots)
}
\arguments{
  \item{lower}{ the vector of lower limits of length n.}
  \item{upper}{ the vector of upper limits of length n.}
  \item{mean}{ the mean vector of length n.}
  \item{corr}{ the correlation matrix of dimension n.}
  \item{sigma}{ the covariance matrix of dimension n less than 1000. Either \code{corr} or
                \code{sigma} can be specified. If \code{sigma} is given, the
                problem is standardized internally. If \code{corr} is given,
                it is assumed that appropriate standardization was performed
                by the user.  If neither \code{corr} nor
                \code{sigma} is given, the identity matrix is used
                for \code{sigma}. }
  \item{algorithm}{ an object of class \code{\link{GenzBretz}},
                    \code{\link{Miwa}} or \code{\link{TVPACK}}
                    specifying both the algorithm to be used as well as
                    the associated hyper parameters.}
  \item{keepAttr}{\code{\link{logical}} indicating if
    \code{\link{attributes}} such as \code{error} and \code{msg} should be
    attached to the return value.  The default, \code{TRUE} is back compatible.}
  \item{seed}{an object specifying if and how the random number generator
          should be initialized, see \code{\link[stats]{simulate}}.}
  \item{\dots}{ additional parameters (currently given to \code{\link{GenzBretz}} for
              backward compatibility issues). }
}
\details{

This program involves the computation of
multivariate normal probabilities with arbitrary correlation matrices.
It involves both the computation of singular and nonsingular
probabilities. The implemented methodology is described in
\bibcitet{mvtnorm::numerical-:1992} and \bibcitet{mvtnorm::comparison:1993} for algorithm
\code{GenzBretz}, in \bibcitet{mvtnorm::Miwa+Hayter+Kuriki:2003} for algorithm
\code{Miwa}, useful up to dimension 20, and \bibcitet{mvtnorm::Genz:2004}
for the \code{TVPACK} algorithm, which covers 2- and 3-dimensional problems
for semi-infinite integration regions.

Note the default algorithm GenzBretz is randomized and hence slightly depends on
\code{\link{.Random.seed}} and that both \code{-Inf} and \code{+Inf} may
be specified in \code{lower} and \code{upper}.  For more details see
\code{\link{pmvt}}.

The multivariate normal
case is treated as a special case of \code{\link{pmvt}} with \code{df=0} and
univariate problems are passed to \code{\link{pnorm}}.

The multivariate normal density and random deviates are available using
\code{\link{dmvnorm}} and \code{\link{rmvnorm}}.

\code{pmvnorm} is based on original implementations by Alan Genz, Frank
Bretz, and Tetsuhisa Miwa developed for computing accurate approximations to
the normal integral. Users interested in computing log-likelihoods involving
such normal probabilities should consider function \code{\link{lpmvnorm}},
which is more flexible and efficient for this task and comes with the
ability to evaluate score functions.

An overview is available from \bibcitet{mvtnorm::Genz_Bretz_2009}.

}
\value{
  The evaluated distribution function is returned, if \code{keepAttr} is true, with attributes
  \item{error}{estimated absolute error}
  \item{msg}{status message(s).}
  \item{algorithm}{a \code{\link{character}} string with \code{class(algorithm)}.}
}
\references{\bibshow{*}}
\seealso{\code{\link{qmvnorm}} for quantiles and \code{\link{lpmvnorm}} for
log-likelihoods.}

\examples{

n <- 5
mean <- rep(0, 5)
lower <- rep(-1, 5)
upper <- rep(3, 5)
corr <- diag(5)
corr[lower.tri(corr)] <- 0.5
corr[upper.tri(corr)] <- 0.5
prob <- pmvnorm(lower, upper, mean, corr)
print(prob)

stopifnot(pmvnorm(lower=-Inf, upper=3, mean=0, sigma=1) == pnorm(3))

a <- pmvnorm(lower=-Inf,upper=c(.3,.5),mean=c(2,4),diag(2))

stopifnot(round(a,16) == round(prod(pnorm(c(.3,.5),c(2,4))),16))

a <- pmvnorm(lower=-Inf,upper=c(.3,.5,1),mean=c(2,4,1),diag(3))

stopifnot(round(a,16) == round(prod(pnorm(c(.3,.5,1),c(2,4,1))),16))

# Example from R News paper (original by Genz, 1992):

m <- 3
sigma <- diag(3)
sigma[2,1] <- 3/5
sigma[3,1] <- 1/3
sigma[3,2] <- 11/15
pmvnorm(lower=rep(-Inf, m), upper=c(1,4,2), mean=rep(0, m), corr=sigma)

# Correlation and Covariance

a <- pmvnorm(lower=-Inf, upper=c(2,2), sigma = diag(2)*2)
b <- pmvnorm(lower=-Inf, upper=c(2,2)/sqrt(2), corr=diag(2))
stopifnot(all.equal(round(a,5) , round(b, 5)))

}
\keyword{distribution}
