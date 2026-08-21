\name{hypergeo_A_nonpos_int}
\alias{hypergeo_A_nonpos_int}
\alias{hypergeo_AorB_nonpos_int}
\title{Hypergeometric functions for integer arguments}
\description{
  Hypergeometric functions for A and/or B being integers
}
\usage{
hypergeo_A_nonpos_int(A, B, C, z, tol = 0)
hypergeo_AorB_nonpos_int(A, B, C, z, tol = 0)
}
\arguments{
  \item{A,B,C}{Parameters for the hypergeometric function}
  \item{z}{Primary complex argument}
  \item{tol}{tolerance}
}
\details{

  The \dQuote{point} of these functions is that if A and C (or B and C)
  are identical nonpositive integers, a warning needs to be given
  because the function is defined as the appropriate limit and one needs
  to be sure that both A and C approach that limit at the same speed.
  
  Function \code{hypergeo_AorB_nonpos_int()} is a convenience wrapper
  for \code{hypergeo_A_nonpos_int()}.
}
\references{
  M. Abramowitz and I. A. Stegun 1965. \emph{Handbook of
    mathematical functions}.  New York: Dover
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{hypergeo}}}
\examples{

jjR1 <- hypergeo(-4, pi, 2.2 , 1+6i)
jjR2 <- hypergeo(pi, -4, 2.2 , 1+6i)  # former bug
jjM <- 3464.1890402837334002-353.94143580568566281i  # value given by Mathematica 


}
\keyword{math}
