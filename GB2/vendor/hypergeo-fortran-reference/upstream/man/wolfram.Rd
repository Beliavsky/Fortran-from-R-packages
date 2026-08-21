\name{wolfram}
\alias{wolfram}
\alias{w07.23.06.0026.01}
\alias{w07.23.06.0026.01_bit1}
\alias{w07.23.06.0026.01_bit2}
\alias{w07.23.06.0026.01_bit3_a}
\alias{w07.23.06.0026.01_bit3_b}
\alias{w07.23.06.0026.01_bit3_c}
\alias{w07.23.06.0029.01}
\alias{w07.23.06.0031.01}
\alias{w07.23.06.0031.01_bit1}
\alias{w07.23.06.0031.01_bit2}
\title{Various functions taken from the Wolfram Functions Site}
\description{
  Various functions taken from the Wolfram Functions Site
}
\usage{
w07.23.06.0026.01(A, n, m, z, tol = 0, maxiter = 2000, method = "a")
w07.23.06.0026.01_bit1(A, n, m, z, tol = 0)
w07.23.06.0026.01_bit2(A, n, m, z, tol = 0, maxiter = 2000)
w07.23.06.0026.01_bit3_a(A, n, m, z, tol = 0)
w07.23.06.0026.01_bit3_b(A, n, m, z, tol = 0)
w07.23.06.0026.01_bit3_c(A, n, m, z, tol = 0)
w07.23.06.0029.01(A, n, m, z, tol = 0, maxiter = 2000)
w07.23.06.0031.01(A, n, m, z, tol = 0, maxiter = 2000)
w07.23.06.0031.01_bit1(A, n, m, z, tol = 0, maxiter = 2000)
w07.23.06.0031.01_bit2(A, n, m, z, tol = 0, maxiter = 2000)
}
\arguments{
  \item{A}{Parameter of hypergeometric function}
  \item{m,n}{Integers}
  \item{z}{Primary complex argument}
  \item{tol,maxiter}{Numerical arguments as per \code{genhypergeo()}}
  \item{method}{Character, specifying method to  be used}
}
\details{

  The \code{method} argument is described at \code{f15.3.10}.  All
  functions' names follow the conventions in
  \code{Hypergeometric2F1.pdf}.

  \itemize{
    \item Function \code{w07.23.06.0026.01(A, n, m, z)} returns
    \eqn{{}_2F_1(A,A+n,A+m,z)}{hypergeo(A,A+n,A+m,z)} where \eqn{m} and
    \eqn{n} are nonnegative integers with \eqn{m\geq n}{m>=n}.
    \item Function \code{w07.23.06.0029.01(A, n, m, z)} returns
    \eqn{{}_2F_1(A,A+n,A-m,z)}{hypergeo(A,A+n,A-m,z)}.
    \item Function \code{w07.23.06.0031.01(A, n, m, z)} returns
    \eqn{{}_2F_1(A,A+n,A+m,z)}{hypergeo(A,A+n,A-m,z)} with \eqn{m\leq
    n}{m<=n}. 
  } 
}
\note{
  These functions use the \code{psigamma()} function which does not yet
  take complex arguments; this means that complex values for \code{A}
  are not supported.  I'm working on it.
  }
\references{
  \code{http://functions.wolfram.com/Hypergeometric2F1.pdf}
}
\author{Robin K. S. Hankin}
\seealso{\code{\link{f15.3.10}},\code{\link{hypergeo}}}
\examples{

# Here we catch some answers from Maple (jjM) and compare it with R's:


jjM <- 0.95437201847068289095 + 0.80868687461954479439i # Maple's answer
jjR <- w07.23.06.0026.01(A=1.1 , n=1 , m=4 , z=1+1i) 
# [In practice, one would type 'hypergeo(1.1, 2.1, 5.1, 1+1i)']

stopifnot(Mod(jjM - jjR) < 1e-10)


jjM <- -0.25955090546083991160e-3 - 0.59642767921444716242e-3i
jjR <- w07.23.06.0029.01(A=4.1 , n=1 , m=1 , z=1+4i)
# [In practice, one would type 'hypergeo(4.1, 3.1, 5.1, 1+1i)']

stopifnot(Mod(jjM - jjR) < 1e-15)

jjM <- 0.33186808222278923715e-1 - 0.40188208572232037363e-1i
jjR <- w07.23.06.0031.01(6.7,2,1,2+1i)
# [In practice, one would type 'hypergeo(6.7, 8.7, 7.7, 2+1i)']
stopifnot(Mod(jjM - jjR) < 1e-10)
}
\keyword{math}
