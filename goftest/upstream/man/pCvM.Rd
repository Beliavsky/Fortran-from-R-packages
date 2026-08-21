\name{pCvM}
\alias{pCvM}
\alias{qCvM}
\title{
  Null Distribution of Cramer-von Mises Test Statistic
}
\description{
  \code{pCvM} computes the cumulative distribution function,
  and \code{qCvM} computes the quantile function,
  of the null distribution of the
  \ifelse{latex}{\out{Cram\'er}}{Cramer}-von Mises test
  statistic.
}
\usage{
  pCvM(q, n = Inf, lower.tail = TRUE)
  qCvM(p, n = Inf, lower.tail = TRUE)
}
\arguments{
  \item{q}{
    Numeric vector of quantiles (values for which the
    cumulative probability is required).
  }
  \item{p}{
    Numeric vector of probabilities.
  }
  \item{n}{
    Integer. Sample size for the
    \ifelse{latex}{\out{Cram\'er}}{Cramer}-von Mises test.
  }
  \item{lower.tail}{
    Logical. If \code{TRUE} (the default),
    probabilities are \eqn{P(X \le q)}{P(X <= q)},
    and otherwise they are \eqn{P(X > q)}.
  }
}
\details{
  For finite \code{n} the cumulative distribution function is
  approximated by the first order expansion
  \eqn{V(x) + \psi_1(x)/n}{V(x) + psi1(x)/n},
  equation (1.8) of
  \ifelse{latex}{\out{Cs\"org\"o}}{Csorgo} and Faraway (1996).

  \code{qCvM} uses \code{\link[stats]{uniroot}} to find the
  quantiles.
}
\value{
  A numeric vector of the same length as \code{p} or \code{q}.
}
\references{
  \ifelse{latex}{\out{Cs\"org\H{o}}}{Csorgo}, S. and Faraway, J.J. (1996)
  The exact and asymptotic distributions of
  \ifelse{latex}{\out{Cram\'er}}{Cramer}-von Mises statistics.
  \emph{Journal of the Royal Statistical Society, Series B}
  \bold{58}, 221--234.
}
\author{
  Original Matlab code by Julian Faraway,
  translated to \R by Adrian Baddeley.
}
\seealso{
 \code{\link{cvm.test}}
}
\examples{
  pCvM(1.1, n=5)
  pCvM(1.1)

  qCvM(0.5, n=5)
  qCvM(0.5)
}
\keyword{distribution}
\keyword{htest}
