\name{bigz_operators}
\alias{add.bigz}
\alias{+.bigz}
\alias{sub.bigz}
\alias{-.bigz}
\alias{mul.bigz}
\alias{*.bigz}
\alias{div.bigz}
\alias{/.bigz}
\alias{mod.bigz}
\alias{\%\%.bigz}
\alias{\%/\%.bigz}
\alias{divq.bigz}
\alias{abs.bigz}
\alias{inv.bigz}
\alias{inv}
\alias{pow.bigz}
\alias{pow}
\alias{^.bigz}
\alias{log.bigz}
\alias{log2.bigz}
\alias{log10.bigz}
\title{Basic Arithmetic Operators for Large Integers ("bigz")}
\description{
  Addition, substraction, multiplication, (integer) division,
  remainder of division, multiplicative inverse,
  power and logarithm functions.
}
\usage{
add.bigz(e1, e2)
sub.bigz(e1, e2 = NULL)
mul.bigz(e1, e2)
div.bigz(e1, e2)
divq.bigz(e1,e2) ## ==  e1 \%/\% e2
mod.bigz(e1, e2) ## ==  e1 \%\%  e2
\method{abs}{bigz}(x)

inv.bigz(a, b,...)## == (1 / a) (modulo b)
% inv(a, ...)
pow.bigz(e1, e2,...)## == e1 ^ e2
% pow(a, ...)
\S3method{log}{bigz}(x, base=exp(1))
\S3method{log2}{bigz}(x)
\S3method{log10}{bigz}(x)
}

\arguments{
  \item{x}{bigz, integer or string from an integer}
  \item{e1, e2, a,b}{bigz, integer or string from an integer}
  \item{base}{base of the logarithm; base e as default}
  \item{...}{Additional parameters}
}
\details{
  Operators can be used directly when objects are of class bigz:
  a + b, log(a), etc.

  For details about the internal modulus state, \emph{and} the rules
  applied for arithmetic operations on big integers with a modulus,
  see the \code{\link{bigz}} help page.

  \code{a / b}  \eqn{=}  \code{div(a,b)} returns a rational number
  unless the operands have a (matching) modulus where
  \code{a * b^-1} results.
  \cr
  \code{a \%/\% b} (or, equivalently, \code{divq(a,b)}) returns the
  quotient of simple \emph{integer} division (with truncation towards zero),
  possibly re-adding a modulus at the end (but \emph{not} using a
  modulus like in \code{a / b}).

  \code{r <- inv.bigz(a, m)}, the multiplicative inverse of
  \code{a} modulo \eqn{m}, corresponds to \code{1/a} or \code{a ^-1}
  from above \emph{when} \code{a} has modulus \code{m}.  Note that
  \eqn{a} not always has an inverse modulo \eqn{m}, in which case
  \code{r} will be \code{\link{NA}} with a warning that can be turned
  off via \preformatted{options("gmp:warnNoInv" = FALSE)}.
}
\value{
  Apart from \code{/} (or \code{div}), where rational numbers
  (\code{\link{bigq}}) may result, these functions return an object of
  class \code{"bigz"}, representing the result of the arithmetic
  operation.
}

\references{The GNU MP Library, see \url{https://gmplib.org}}

\author{Immanuel Scholz and Antoine Lucas}

\examples{
# 1+1=2
as.bigz(1) + 1
as.bigz(2)^10
as.bigz(2)^200

# if my.large.num.string is set to a number, this returns the least byte
(my.large.num.string <- paste(sample(0:9, 200, replace=TRUE), collapse=""))
mod.bigz(as.bigz(my.large.num.string), "0xff")

# power exponents can be up to MAX_INT in size, or unlimited if a
# bigz's modulus is set.
pow.bigz(10,10000)

## Modulo 11,   7 and 8 are inverses :
as.bigz(7, mod = 11) * 8 ## ==>  1  (mod 11)
inv.bigz(7, 11)## hence, 8
a <- 1:10
(i.a <- inv.bigz(a, 11))
d <- as.bigz(7)
a \%/\% d  # = divq(a, d)
a \%\%  d  # = mod.bigz (a, d)

(ii <- inv.bigz(1:10, 16))
## with 5 warnings (one for each NA)
op <- options("gmp:warnNoInv" = FALSE)
i2 <- inv.bigz(1:10, 16) # no warnings
(i3 <- 1 / as.bigz(1:10, 16))
i4 <- as.bigz(1:10, 16) ^ -1
stopifnot(identical(ii, i2),
	  identical(as.bigz(i2, 16), i3),
	  identical(i3, i4))
options(op)# revert previous options' settings

stopifnot(inv.bigz(7, 11) == 8,
          all(as.bigz(i.a, 11) * a == 1),
          identical(a \%/\% d, divq.bigz(1:10, 7)),
          identical(a \%\%  d, mod.bigz (a, d))
 )

}

\keyword{arith}
%\keyword{bigz}
