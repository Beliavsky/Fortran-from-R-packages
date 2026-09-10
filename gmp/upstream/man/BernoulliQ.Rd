\name{BernoulliQ}
\alias{BernoulliQ}
\title{Exact Bernoulli Numbers}
\description{%% ../R/Stirling-n-etc.R :
  Return the \eqn{n}-th Bernoulli number \eqn{B_n}, (or \eqn{B_n^+}{Bn+},
  see the reference), where \eqn{B_1 = + \frac 1 2}{B_1 = + 1/2}.
}
\usage{
BernoulliQ(n, verbose = getOption("verbose", FALSE))
}
\arguments{
  \item{n}{integer \emph{vector}, \eqn{n \ge 0}{n >= 0}.}
  \item{verbose}{logical indicating if computation should be traced.}
}
\value{
  a big rational (class \code{\link[=bigq-class]{"bigq"}}) vector of the
  Bernoulli numbers \eqn{B_n}.
}
\references{\url{https://en.wikipedia.org/wiki/Bernoulli_number}
}
\author{Martin Maechler}
\seealso{
  \code{\link[Rmpfr]{Bernoulli}} in \CRANpkg{Rmpfr} in arbitrary precision
  via Riemann's \eqn{\zeta}{zeta} function.
  % \code{\link[DPQ]{Bern}(n)}
  \code{Bern(n)} in \CRANpkg{DPQ} uses standard (double precision)
  \R arithmetic for the n-th Bernoulli number.
}
\examples{
(Bn0.10 <- BernoulliQ(0:10))
}
