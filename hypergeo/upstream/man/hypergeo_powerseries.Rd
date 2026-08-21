\name{hypergeo_powerseries}
\alias{hypergeo_powerseries}
\alias{hypergeo_taylor}
\alias{hypergeo_general}
\title{The hypergeometric function as determined by power series}
\description{
  The hypergeometric function as determined by infinite
  (\code{hypergeo_powerseries()}) or finite (\code{hypergeo_taylor()})
  power series
}
\usage{
hypergeo_powerseries(A, B, C, z, tol = 0, maxiter = 2000)
}
\arguments{
  \item{A,B,C}{Parameters of the hypergeometric function}
  \item{z}{Primary complex argument}
  \item{tol,maxiter}{Numerical arguments}
}
\details{
  Function \code{hypergeo_powerseries()} is the primary decision-making
  function of the package.  It is this function that detects degenerate
  cases of the three parameters and dispatches accordingly.
  Non-degenerate cases are sent to function \code{hypergeo_general()}.

  Function \code{hypergeo_taylor()} deals with cases where the
  hypergeometric function is a polynomial.
}
\references{ 
  M. Abramowitz and I. A. Stegun 1965. \emph{Handbook of
    mathematical functions}.  New York: Dover
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{hypergeo}},\code{\link{genhypergeo}}}
\examples{


jjR <- hypergeo(pi,-4,2.2,1+5i)
jjM <- 1670.8287595795885335 - 204.81995157365381258i


}
\keyword{math}
