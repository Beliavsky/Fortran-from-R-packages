\name{shanks}
\alias{shanks}
\alias{hypergeo_shanks}
\alias{genhypergeo_shanks}
\title{Evaluation of the hypergeometric function using Shanks's method}
\description{
Evaluation of the hypergeometric function using Shanks transformation of
successive sums
}
\usage{
hypergeo_shanks(A,B,C,z,maxiter=20)
genhypergeo_shanks(U,L,z,maxiter=20)
shanks(Last,This,Next)
}
\arguments{
  \item{A,B,C}{Parameters (real or complex)}
  \item{U,L}{Upper and lower vectors}
  \item{z}{Primary complex argument}
  \item{maxiter}{Maximum number of iterations}
  \item{Last,This,Next}{Three successive convergents}
}
\details{

  The Shanks transformation of successive partial sums is

  \deqn{S(n)=\frac{A_{n+1}A_{n-1}-A_n^2}{A_{n+1}-2A_n+A_{n-1}}}


  and if the \eqn{A_n} tend to a limit then the sequence \eqn{S(n)}
  often converges more rapidly than \eqn{A_n}.  However, the denominator
  is susceptible to catastrophic rounding under fixed-precision
  arithmetic and it is difficult to know when to stop iterating.

}
\references{ 
  \itemize{
    \item Shanks, D. (1955). \dQuote{Non-linear transformation of
      divergent and slowly convergent sequences}, \emph{Journal of
      Mathematics and Physics} 34:1-42
  }
}
\author{Robin K. S. Hankin}
\note{
The 
}
\seealso{\code{\link{buhring}}}
\examples{

hypergeo_shanks(1/2,1/3,pi,z= 0.1+0.2i)
}
\keyword{math}

