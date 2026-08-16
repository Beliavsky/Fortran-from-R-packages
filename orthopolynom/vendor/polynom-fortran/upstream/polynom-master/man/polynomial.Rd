\name{polynomial}
\alias{polynomial}
\alias{as.polynomial}
\alias{is.polynomial}
\alias{as.character.polynomial}
\alias{print.polynomial}
\alias{coef.polynomial}
\title{Polynomials}
\description{
  Construct, coerce to, test for, and print polynomial objects.
}
\usage{
polynomial(coef = c(0, 1))
as.polynomial(p)
is.polynomial(p)

\method{as.character}{polynomial}(x, decreasing = FALSE, \dots)
\method{print}{polynomial}(x, digits = getOption("digits"), decreasing = FALSE, \dots)
}
\arguments{
  \item{coef}{numeric vector, giving the polynomial coefficients in
    \emph{in}creasing order.}
  \item{p}{an arbitrary \R object.}
  \item{x}{a \code{polynomial} object.}
  \item{decreasing}{a logical specifying the order of the terms;
    in increasing (default) or decreasing powers.}
  \item{digits}{the number of significant digits to use for printing.}
  \item{\dots}{potentially further arguments passed to and from other methods.}
}
\details{
  \code{polynomial} constructs a polynomial from its coefficients,
  i.e., \code{p[1:k]} specifies the polynomial
  \deqn{p_1 + p_2 x + p_3 x^2 + \dots + p_k x^{k-1}.}{p[1] + p[2]* x +
    p[3]* x^2 + ... + p[k]* x^(k-1).}
  Internally, polynomials are simply numeric coefficient vectors of
  class \code{"polynomial"}.  Several useful methods are available for
  this class, such as coercion to character (\code{as.character()}) and
  function (\code{\link{as.function.polynomial}}), extraction of
  the coefficients (\code{coef()}), printing (using \code{as.character}),
  plotting (\code{\link{plot.polynomial}}), and computing sums and
  products of arbitrarily many polynomials.

  \code{as.polynomial} tries to coerce its arguments to a polynomial.

  \code{is.polynomial} tests whether its argument is a polynomial (in
  the sense that it has class \code{"polynomial"}.
}
\examples{
polynomial(1:4)
p <- as.polynomial(c(1,0,3,0))
p
print(p, decreasing = TRUE)
stopifnot(coef(p) == c(1,0,3))

polynomial(c(2,rep(0,10),1))
}
\keyword{symbolmath}
