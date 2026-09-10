\name{bigz}
\title{Large Sized Integer Values}
\alias{bigz}
\alias{bigz-class}
\alias{as.bigz}
\alias{is.bigz}
\alias{as.character.bigz}
\alias{as.double.bigz}
\alias{print.bigz}
\alias{is.na.bigz}
\alias{NA_bigz_}
\alias{c_bigz}
% these two are C symbols, will be removed eventually:
\alias{biginteger_as}
\alias{biginteger_as_character}
\description{
  Class \code{"bigz"} encodes arbitrarily large integers (via GMP).
  A simple S3 class (internally a \code{\link{raw}} vector), it has been
  registered as formal (S4) class (via \code{\link{setOldClass}}), too.
}
\usage{
as.bigz(a, mod = NA)
NA_bigz_
\S3method{as.character}{bigz}(x, b = 10, \dots)
% \S3method{as.double}{bigz}(x,\dots)
is.bigz(x)
\S3method{is.na}{bigz}(x)
\S3method{print}{bigz}(x, quote=FALSE, initLine = is.null(modulus(x)), \dots)
c_bigz(L)
}
\arguments{
  \item{a}{either \code{\link{integer}}, \code{\link{numeric}} (i.e.,
    \code{\link{double}}) or \code{\link{character}} vector.

    If character: the strings either start with \code{0x} for hexadecimal, \code{0b} for
    binary, \code{0} for octal, or without a \code{0*} prefix for decimal values.
    Formatting errors are signalled as with \code{\link{stop}}.}
  \item{b}{base: from 2 to 36}
  \item{x}{a \dQuote{big integer number} (vector), of class \code{"bigz"}.}
  \item{\dots}{additional arguments passed to methods}
  \item{mod}{an integer, numeric, string or bigz of the internal modulus, see below.}
  \item{quote}{(for printing:) logical indicating if the numbers
    should be quoted (as characters are);  the default used to be
    \code{TRUE} (implicitly) till 2011.}
  \item{initLine}{(for printing:) logical indicating if an
    \bold{init}ial line (with the class and length or dimension) should
    be printed.  The default prints it for those cases where the class
    is not easily discernable from the print output.}
  \item{L}{a \code{\link{list}} where each element contains \code{"bigz"}
    numbers, for \code{c_bigz()}, this allows something like an
    \code{\link{sapply}()} for \code{"bigz"} vectors, see \code{sapplyZ()}
    in the examples.}
}
\value{
  An \R object of (S3) class \code{"bigz"}, representing the argument
  (\code{x} or \code{a}).
}
\details{
  Bigz's are integers of arbitrary, but given length (means: only
  restricted by the host memory).  Basic arithmetic operations can be performed
  on bigzs as addition, subtraction, multiplication, division, modulation
  (remainder of division), power, multiplicative inverse, calculating of the
  greatest common divisor, test whether the integer is prime and other operations
  needed when performing standard cryptographic operations.

  For a review of basic arithmetics, see \code{\link{add.bigz}}.

  Comparison are supported, i.e., \code{"=="}, \code{"!="},
  \code{"<"}, \code{"<="}, \code{">"}, and \code{">="}.

  \code{NA_bigz_} is computed on package load time as \code{as.bigz(NA)}.

  Objects of class \code{"bigz"} may have a \dQuote{modulus}, accessible
  via \code{\link{modulus}()}, currently as an attribute \code{mod}.
  When the object has such a modulus \eqn{m}, arithmetic is performed
  \emph{\dQuote{modulo m}}, mathematically \dQuote{within the
    ring \eqn{Z/mZ}}.  For many operations, this means
\preformatted{   result <- mod.bigz(result, m)  ## == result \%\% m
}
  is called after performing the arithmetic operation and the result
  will have the attribute \code{mod} set accordingly.
  This however does not apply, e.g., for \code{/}, where
  \eqn{a / b := a b^{-1}}{a / b =  a b^(-1)} and
  \eqn{b^{-1}}{b^(-1)} is the \emph{multiplicate inverse} of \eqn{b}
  with respect to ring arithmetic, or \code{\link{NA}} with a warning
  when the inverse does not exist.  The warning can be turned off via
  \code{options("gmp:warnModMismatch" = FALSE)}% FIXME: show+check option()

  Powers of bigzs can only be performed, if either a modulus is going to
  be applied to the result bigz or if the exponent fits into an integer
  value.  So, if you want to calculate a power in a finite group
  (\dQuote{modulo c}), for large \eqn{c} do not use
  \code{a ^ b \%\% c}, but rather \code{as.bigz(a,c) ^ b}.

  The following rules for the result's modulus apply when performing
  arithmetic operations on \code{bigz}s:
  \itemize{
    \item If none of the operand has a modulus set, the result will not have a modulus.
    \item If both operands have a different modulus, the result will not have a
    modulus, except in case of \code{\link{mod.bigz}}, where the second operand's
    value is used.
    \item If only one of the operands has a modulus or both have a common (the
    same), it is set and used for the arithmetic operations, except in
    case of \code{mod.bigz}, where the second operand's value is used.
  }
}

\note{
  \preformatted{
    x <- as.bigz(1234567890123456789012345678901234567890)
  }
  will not work as \R converts the number to a double, losing precision
  and only then convert to a \code{"bigz"} object.

  Instead, use the syntax
  \preformatted{
    x <- as.bigz("1234567890123456789012345678901234567890")
  }
}

\references{The GNU MP Library, see \url{https://gmplib.org}}
\author{Immanuel Scholz}

\examples{
## 1+1=2
a <- as.bigz(1)
a + a

## Two non-small Mersenne primes:
two <- as.bigz(2)
p1 <- two^107 -1 ; isprime(p1); p1
p2 <- two^127 -1 ; isprime(p2); p2

stopifnot( is.na(NA_bigz_) )

## Calculate c = x^e mod n
## --------------------------------------------------------------------
x <- as.bigz("0x123456789abcdef") # my secret message
e <- as.bigz(3) # something smelling like a dangerous public RSA exponent
(n <- p1 * p2) #  a product of two primes
as.character(n, b=16)# as both primes were Mersenne's..

## recreate the three numbers above [for demo below]:
n. <- n; x. <- x; e. <- e # save
Rev <- function() { n <<- n.; x <<- x.; e <<- e.}

# first way to do it right
modulus(x) <- n
c <- x ^ e ; c ; Rev()

# similar second way (makes more sense if you reuse e) to do it right
modulus(e) <- n
c2 <- x ^ e
stopifnot(identical(c2, c), is.bigz(c2)) ; Rev()

# third way to do it right
c3 <- x ^ as.bigz(e, n) ; stopifnot(identical(c3, c))

# fourth way to do it right
c4 <- as.bigz(x, n) ^ e ; stopifnot(identical(c4, c))

# WRONG! (although very beautiful. Ok only for very small 'e' as here)
cc <- x ^ e \%\% n
cc == c

# Return result in hexa
as.character(c, b=16)

# Depict the "S4-class" bigz, i.e., the formal (S4) methods:
if(require("Rmpfr")) # mostly interesting there
  showMethods(class="bigz")

# an  sapply() version that works for big integers "bigz":
sapplyZ <- function(X, FUN, ...) c_bigz(lapply(X, FUN, ...))

# dummy example showing it works (here):
zz <- as.bigz(3)^(1000+ 1:999)
z1 <- sapplyZ(zz, function(z) z^2)
stopifnot( identical(z1, zz^2) )
}
\keyword{arith}
