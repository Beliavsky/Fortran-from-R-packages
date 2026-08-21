\name{quad.form}
\alias{quad.form}
\alias{quadform}
\alias{quad.form.inv}
\alias{quad.form.chol}
\alias{quad3.form}
\alias{quad3.form_ab}
\alias{quad3.form_bc}
\alias{quad3.form.inv}
\alias{quad.tform}
\alias{quad3.tform}
\alias{quad3.tform_ab}
\alias{quad3.tform_bc}
\alias{quad.tform.inv}
\alias{quad.diag}
\alias{quad.tdiag}
\alias{quad3.diag}
\alias{quad3.tdiag}
\alias{quad.trace}
\alias{quad.ttrace}
\alias{cprod}
\alias{tcprod}
\alias{ht}
\alias{cp}
\alias{tcp}
\alias{qf}
\alias{qfi}
\alias{q3}
\alias{q3i}
\alias{q3t}
\alias{qt}
\alias{qti}
\alias{qd}
\alias{qtd}
\alias{q3d}
\alias{q3td}
\alias{qtr}
\alias{qttr}

\title{Evaluate a quadratic form efficiently}
\description{
  
  Given a square matrix \eqn{M} of size \eqn{n\times n}{n*n}, and a
  matrix \eqn{x} of size \eqn{n\times p}{n*p} (or a vector of length
  \eqn{n}), evaluate various quadratic forms.

  The archetype is \code{quad.form(M,x)} for real or complex square
  matrix \code{M} and vector or matrix \code{x}.  This evaluates
  \verb{Conj(t(x)) \%*\% M \%*\% x} but using
  \code{crossprod(crossprod(M,Conj(x)),x)} thus avoiding taking a
  needless transpose.


}
\usage{
quad.form(M, x)
quad.form.inv(M, x)
quad.form.chol(chol, x)
quad.tform(M, x)
quad3.form(M, left, right)
quad3.tform(M, left, right)
quad.tform.inv(M, x)
quad.diag(M, x)
quad.tdiag(M, x)
quad3.diag(M, left, right)
quad3.tdiag(M, left, right)
cprod(x, y)
tcprod(x, y)
ht(x)
}
\arguments{
  \item{M}{Square matrix of size \eqn{n\times n}{n*n}}
  \item{x, y}{Matrix of size \eqn{n\times p}{n*p}, or vector of length \eqn{n}}
  \item{chol}{Lower triangular Cholesky decomposition
    of the quadratic form, see details}
  \item{left, right}{In function \code{quad3.form()}, matrices with
    \eqn{n} rows and arbitrary number of columns}
}
\details{

  \tabular{llll}{
  \code{ht(x)}\tab\eqn{\mathbf{x}^{*}=\overline{\mathbf{x}^{\scriptscriptstyle T}}}{omitted}
  \tab\code{Conj(t(x))}\tab\code{ht(x)}\cr
  \code{cprod(x,y)}\tab\eqn{\mathbf{x}^{*}\mathbf{y}}{omitted}
  \tab\verb{ht(x) \%*\% y}\tab\code{cp()}\cr
  \code{tcprod(x,y)}\tab\eqn{\mathbf{x}\mathbf{y}^{*}}{omitted}
  \tab\verb{x \%*\% ht(y)}\tab\code{tcp()}\cr
  \code{quad.form(M,x)}\tab\eqn{\mathbf{x}^{*}M\mathbf{x}}{omitted}
  \tab\verb{ht(x) \%*\% M \%*\% x}\tab\code{qf()}\cr
  \code{quad.form.inv(M,x)}\tab\eqn{\mathbf{x}^{*}M^{-1}\mathbf{x}}{omitted}
  \tab\verb{ht(x) \%*\% solve(M) \%*\% x}\tab\code{qfi()}\cr
  \code{quad.tform(M,x)}\tab\eqn{\mathbf{x}M\mathbf{x}^{*}}{omitted}
  \tab\verb{x \%*\% A \%*\% ht(x)}\tab\code{qt()}\cr
  \code{quad.tform.inv(M,x)}\tab\eqn{\mathbf{x}M^{-1}\mathbf{x}^{*}}{omitted}
  \tab\verb{x \%*\% solve(M) \%*\% ht(x)}\tab\code{qti()}\cr
  \code{quad3.form(M,l,r)}\tab\eqn{\mathbf{l}^{*}M\mathbf{r}}{omitted}
  \tab\verb{t(l) \%*\% M \%*\% r}\tab\code{q3()}\cr
  \code{quad3.form.inv(M,l,r)}\tab\eqn{\mathbf{l}^{*}M^{-1}\mathbf{r}}{omitted}
  \tab\verb{t(l) \%*\% solve(M) \%*\% r}\tab\code{q3i()}\cr
  \code{quad3.tform(M,l,r)}\tab\eqn{\mathbf{l}M\mathbf{r}^{*}}{omitted}
  \tab\verb{l \%*\% M \%*\% t(r)}\tab\code{q3t()}\cr
  \code{quad.diag(M,x)}\tab\eqn{\operatorname{diag}(\mathbf{x}^{*}M\mathbf{x})}{omitted}
  \tab\verb{diag(quad.form(M,x))}\tab\code{qd()}\cr
  \code{quad.tdiag(M,x)}\tab\eqn{\operatorname{diag}(\mathbf{x}M\mathbf{x}^{*})}{omitted}
  \tab\verb{diag(quad.tform(M,x))}\tab\code{qtd()}\cr
  \code{quad3.diag(M,l,r)}\tab\eqn{\operatorname{diag}(\mathbf{l}^{*}M\mathbf{r})}{omitted}
  \tab\code{diag(quad3.form(M,l,r))}\tab\code{q3d()}\cr
  \code{quad3.tdiag(M,l,r)}\tab\eqn{\operatorname{diag}(\mathbf{l}M\mathbf{r}^{*})}{omitted}
  \tab\code{diag(quad3.tform(M,l,r))}\tab\code{q3td()}\cr
  \code{quad.trace(M,x)}\tab\eqn{\operatorname{tr}(\mathbf{x}^{*}M\mathbf{x})}{omitted}
  \tab\code{tr(quad.form(M,x))}\tab\code{qt()}\cr
  \code{quad.ttrace(M,x)}\tab\eqn{\operatorname{tr}(\mathbf{x}M\mathbf{x}^{*})}{omitted}
  \tab\code{tr(quad.tform(M,x))}\tab\code{qtt()}
}

\strong{In the above, \eqn{\mathbf{x}^{*}}{ht(x)} denotes the
  \emph{complex conjugate} of the transpose, also known as the Hermitian
  transpose} (this only matters when considering complex numbers).

  These various functions generally avoid taking needless expensive
  transposes in favour of using nested \code{crossprod()} and
  \code{tcrossprod()} calls.  For example, the \dQuote{meat} of
  \code{quad.form()} is just \code{crossprod(crossprod(M,Conj(x)),x)}.

  Functions such as \code{quad.form.inv()} avoid taking a matrix
  inverse.  The meat of \code{quad.form.inv()}, for example, is
  \code{cprod(x, solve(M, x))}.  Many people have stated things like
  \dQuote{Never invert a matrix unless absolutely necessary}.  But I
  have \emph{never} seen a case where \code{quad.form.inv(M,x)} is
  faster than \code{quad.form(solve(M),x)}.

  One motivation for the package is to return consistent results with
  complex arguments.  Note, for example, that
  \code{base::crossprod(x,y)} evaluates \verb{t(x) \%*\% y} and not, as
  one would almost always want, \verb{Conj(t(x)) \%*\% y}.  Function
  \code{cprod()}, unlike \code{crossprod()}, is consistent and returns
  \verb{Conj(t(x)) \%*\% y} [or \verb{ht(x) \%*\% y}]; internally it is
  essentially \verb{crossprod(Conj(x), y)}.
  
  Function \code{quad.form.chol()} interprets argument \code{chol} as
  the lower triangular Cholesky decomposition of the quadratic form.
  Remember that \code{M.lower \%*\% M.upper == M}, and \code{chol()}
  returns the upper triangular matrix, so one needs to use the transpose
  \code{t(chol(M))} in calls.  If the Cholesky decomposition of \code{M}
  is available, then using \code{quad.form.chol()} and supplying
  \code{chol(M)} should generally be faster (for large matrices) than
  calling \code{quad.form()} and using \code{M} directly.  The time
  saving is negligible for matrices smaller than about \eqn{50\times
  50}{50*50}, even if the overhead of computing the decomposition is
  ignored.

  Functions \code{quad3.foo()} take three arguments: a matrix \code{M}
  and two other vectors \code{l} and \code{r} [or \code{left} and
  \code{right}].  For these functions, \code{M} is not necessarily
  square although of course the matrices have to be compatible.

  Functions \code{quad3.form_ab()} and \code{quad3.form_bc()} are helper
  functions not really intended for the end-user.  They return
  mathematically identical results but differ in the bracketing order of
  their operations: \code{quad3.form_ab(M,l,r)} returns
  \eqn{\left(\mathbf{l}^*M\right)\mathbf{r}}{omitted} and
  \code{quad3.form_bc(M,l,r)} returns
  \eqn{\mathbf{l}^{*}\left(M\mathbf{r}\right)}{omitted}.  The mnemonic
  for their names is derived from the first multiplication when
  calculating \eqn{(ab)c} and \eqn{a(bc)}.  Note that
  \code{quad3.form_ab(M,l,r)} returns
  \verb{crossprod(crossprod(M,Conj(l)),r)} rather than the
  mathematically equivalent \verb{cprod(cprod(M,l),r)} on efficiency
  grounds (only a single conjugate is taken).

  Function \code{quad3.form()} dispatches to either
  \code{quad3.form_ab()} or \code{quad3.form_bc()} depending on the
  dimensions of its argument as per the efficiency discussion at
  \code{inst/quadform3test.Rmd}.  Similar considerations apply to
  \code{quad3.tform()}, \code{quad3.tform_ab()}, and
  \code{quad3.tform_bc()}.

  Terse forms [\code{qf()} for \code{quad.form()}, \code{qti()} for
  \code{quad.tform.inv()}, etc] are provided for the \proglang{perl}
  golfers among us.

}

\value{
  Generally, return a (dropped) matrix, real or complex as appropriate
}

\author{Robin K. S. Hankin}

\note{

 These functions are used extensively in the \CRANpkg{emulator} and
 \CRANpkg{calibrator} packages, primarily in the interests of elegant
 code, but also speed.  For the problems I usually consider, the speedup
 (of \code{quad.form(M, x)} over \verb{t(x) \%*\% M \%*\% x}, say) is
 marginal at best.
  
}

\examples{
jj <- matrix(rnorm(80), 20, 4)
M <- crossprod(jj, jj)
M.lower <- t(chol(M))
x <- matrix(rnorm(8), 4, 2)

jj.1 <- t(x) \%*\% M \%*\% x
jj.2 <- quad.form(M, x)
jj.3 <- quad.form.chol(M.lower, x)
print(jj.1)
print(jj.2)
print(jj.3)


## Make two Hermitian positive-definite matrices:
L <- matrix(c(1, 0.1i,-0.1i,1), 2, 2)
LL <- diag(11)
LL[2,1] <- -(LL[1,2] <- 0.1i)

z <- matrix(rnorm(22) + 1i*rnorm(22), 2, 11)

quad.diag(L, z)     # elements real because L is HPD
quad.tdiag(LL, z)   # ditto


## Now consider accuracy:
quad.form(solve(M), x) - quad.form.inv(M, x)  # should be zero
quad.form(M, x) - quad.tform(M, t(x))         # should be zero
quad.diag(M, x) - diag(quad.form(M, x))       # should be zero
diag(quad.form(L, z))   - quad.diag(L, z)     # should be zero
diag(quad.tform(LL, z)) - quad.tdiag(LL, z)   # should be zero
}
\keyword{array}
