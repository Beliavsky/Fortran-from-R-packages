\name{gosper}
\alias{hypergeo_gosper}
\title{Evaluation of the hypergeometric function using Gosper's method}
\description{
  Evaluation of the hypergeometric function using Gosper's method
}
\usage{
hypergeo_gosper(A, B, C, z, tol = 0, maxiter = 2000)
}
\arguments{
  \item{A,B,C}{Parameters (real or complex)}
  \item{z}{Complex argument}
  \item{tol}{tolerance (passed to \code{GCF()})}
  \item{maxiter}{maximum number of iterations}
}
\details{
Gosper provides a three-term recurrence which converges when \eqn{z} is
close to a critical point.

Bill Gosper asserts that the recursion holds for values of \eqn{z} which
are inside the cardioid \code{(sqrt(8)*cos(t)-2*cos(2t),
sqrt(8)*sin(t)-2*sin(2t))} (see examples section).

It is suggested that the recursion should only be used when the
auxiliary parameters A, B,C are all \eqn{\le 12}{<=12} in absolute
value.
}

\references{
  Original email was archived at
  \code{https://www.ma.utexas.edu/pipermail/maxima/2006/000126.html} but
  does not appear there now; and the wayback machine doesn't find it
  either.
}

\author{
  R code by Robin K. S. Hankin, transcribed from maxima code posted by
  Richard Fateman, who credited Bill Gosper
}

\seealso{\code{\link{hypergeo_contfrac}}}
\examples{


hypergeo_gosper(1.1,5.1,3.1,crit())

# Compare MMA: -0.192225 + 0.692328 I

t <- seq(from=0,to=2i*pi,len=100)
plot(exp(t)*(sqrt(8)-exp(t)),asp=1,type='l')
points(crit())


}
\keyword{math}
