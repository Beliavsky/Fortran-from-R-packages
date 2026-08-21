\name{f15.3.10}
\alias{f15.1.1}
\alias{f15.3.10}
\alias{f15.3.10_a}
\alias{f15.3.10_b}
\alias{f15.3.11}
\alias{f15.3.11_bit1}
\alias{f15.3.11_bit2_a}
\alias{f15.3.11_bit2_b}
\alias{f15.3.12}
\alias{f15.3.12_bit1}
\alias{f15.3.12_bit2_a}
\alias{f15.3.12_bit2_b}
\alias{f15.3.13}
\alias{f15.3.13_14}
\alias{f15.3.13_a}
\alias{f15.3.13_b}
\alias{f15.3.14}
\alias{f15.3.14_bit1_a}
\alias{f15.3.14_bit1_b}
\alias{f15.3.14_bit2}
\alias{f15.3.10_11_12}
\title{Transformations of the hypergeometric function}
\description{
Transformations of the hypergeometric function detailed in AMS-55, page
559-560.
}
\usage{
f15.3.10       (A, B,    z, tol = 0, maxiter = 2000, method = "a")
f15.3.10_a     (A, B,    z, tol = 0, maxiter = 2000              )
f15.3.10_b     (A, B,    z, tol = 0, maxiter = 2000              )
f15.3.11       (A, B, m, z, tol = 0, maxiter = 2000, method = "a")
f15.3.11_bit1  (A, B, m, z, tol = 0                              )
f15.3.11_bit2_a(A, B, m, z, tol = 0, maxiter = 2000              )
f15.3.11_bit2_b(A, B, m, z, tol = 0, maxiter = 2000              )
f15.3.12       (A, B, m, z, tol = 0, maxiter = 2000, method = "a")
f15.3.12_bit1  (A, B, m, z, tol = 0                              )
f15.3.12_bit2_a(A, B, m, z, tol = 0, maxiter = 2000              )
f15.3.12_bit2_b(A, B, m, z, tol = 0, maxiter = 2000              )
f15.3.13       (A, C,    z, tol = 0, maxiter = 2000, method = "a")
f15.3.13_a     (A, C,    z, tol = 0, maxiter = 2000              )
f15.3.13_b     (A, C,    z, tol = 0, maxiter = 2000              )
f15.3.14       (A, C, m, z, tol = 0, maxiter = 2000, method = "a")
f15.3.14_bit1_a(A, C, m, z, tol = 0, maxiter = 2000              )
f15.3.14_bit1_b(A, C, m, z, tol = 0, maxiter = 2000              )
f15.3.14_bit2  (A, C, m, z, tol = 0                              )
f15.3.13_14    (A, C, m, z, tol = 0, maxiter = 2000, method = "a")
f15.3.10_11_12 (A, B, m, z, tol = 0, maxiter = 2000, method = "a")
f15.1.1        (A, B, C, z, tol = 0, maxiter = 2000              )
}
\arguments{
  \item{A,B,C}{Parameters of the hypergeometric function}
  \item{m}{Integer linking \code{A}, \code{B}, \code{C} as set out in
AMS-55, page 559,560}
  \item{z}{primary complex argument}
  \item{tol,maxiter}{numerical parameters}
  \item{method}{Length 1 character vector specifying the method.  See details} 
}
\details{

Naming scheme (functions and arguments) follows AMS-55, pages 559-560.

The \code{method} argument to (eg) \code{f15.3.14()} specifies whether
to use \code{psigamma()} directly (method \dQuote{\code{a}}), or the
recurrence 6.3.5 (method \dQuote{\code{b}}).  Press et al recommend
method \dQuote{\code{b}}, presumably on the grounds of execution speed.
I'm not so sure (method \dQuote{\code{a}} seems to be more accurate in
the sense that it returns values closer to those of Maple).

Method \dQuote{\code{c}} means to use a totally dull, slow, direct (but
clearly correct) summation, for the purposes of debugging.  This is
only used for the functions documented under \code{wolfram.Rd}

Functions \code{f15.3.13_14()} and \code{f15.3.10_11_12()} are
convenience wrappers.  For example, function \code{f15.3.13_14()}
dispatches to either \code{f15.3.13()} or \code{f15.3.14()} depending on
the value of \code{m}.
}
\references{
  M. Abramowitz and I. A. Stegun 1965. \emph{Handbook of
    mathematical functions}.  New York: Dover
}
\author{Robin K. S. Hankin}
\note{
  These functions are not really designed to be called by the user: use
  \code{hypergeo()} instead, or \code{hypergeo_cover[123]()} for
  specific cases.

}
\seealso{\code{\link{hypergeo}},\code{\link{wolfram}},\code{\link{hypergeo_cover1}}}
\examples{

f15.3.10_11_12(A=1.1, B=pi, m= +3, z=.1+.1i)
f15.3.10_11_12(A=1.1, B=pi, m= -3, z=.1+.1i)

}
\keyword{math}
