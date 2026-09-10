\name{bigq}
\title{Large sized rationals}
\alias{bigq}
\alias{bigq-class}
\alias{as.bigq}
\alias{is.bigq}
\alias{as.character.bigq}
\alias{as.double.bigq}
\alias{print.bigq}
\alias{is.na.bigq}
\alias{NA_bigq_}
\alias{as.bigz.bigq}
\alias{c_bigq}
\alias{denominator}
\alias{denominator<-}
\alias{numerator}
\alias{numerator<-}
\description{
  Class \code{"bigq"} encodes rationals encoded as ratios of arbitrary
  large integers (via GMP).
  A simple S3 class (internally a \code{\link{raw}} vector), it has been
  registered as formal (S4) class (via \code{\link{setOldClass}}), too.
}
\usage{
as.bigq(n, d = 1)
\S3method{as.character}{bigq}(x, b=10,...)
\S3method{as.double}{bigq}(x,...)
as.bigz.bigq(a, mod=NA)
is.bigq(x)
\S3method{is.na}{bigq}(x)
\S3method{print}{bigq}(x, quote=FALSE, initLine = TRUE, ...)
denominator(x)
numerator(x)
NA_bigq_
c_bigq(L)
}

\arguments{
  \item{n,d}{either integer, numeric or string value
    (String value: either starting with \code{0x} for hexadecimal, \code{0b} for
    binary or without prefix for decimal values.
    Any format error results in \code{0}).
    \code{n} stands for numerator, \code{d} for denominator.}
  \item{a}{an element of class \code{"bigq"}}
  \item{mod}{optional modulus to convert into biginteger}
  \item{x}{a \dQuote{rational number} (vector), of class \code{"bigq"}.}
  \item{b}{base: from 2 to 36}
  \item{...}{additional arguments passed to methods}
  \item{quote}{(for printing:) logical indicating if the numbers
    should be quoted (as characters are);  the default used to be
    \code{TRUE} (implicitly) till 2011.}
  \item{initLine}{(for printing:) logical indicating if an
    \bold{init}ial line (with the class and length or dimension) should
    be printed. % not yet smart default .. The default prints it for
    % those cases where the class is not easily discernable from the
    % print output.
  }
  \item{L}{a \code{\link{list}} where each element contains \code{"bigq"}
    numbers, for \code{c_bigq()}, this allows something like an
    \code{\link{sapply}()} for \code{"bigq"} vectors, see \code{sapplyQ()}
    in the examples below.}
}
\value{
  An \R object of (S3) class \code{"bigq"} representing the parameter value.
}
\details{
  \code{as.bigq(x)} when \code{x} is \code{\link{numeric}} (aka
  \code{\link{double}} precision) calls the \file{GMP} function
  \code{mpq_set_d()} which is documented to be \emph{exact} (every finite double
  precision number is a rational number).

  \code{as.bigz.bigq()} returns the smallest integers not less than the
  corresponding rationals bigq.

  \code{NA_bigq_} is computed on package load time as \code{as.bigq(NA)}.
}

\author{Antoine Lucas}

\examples{
x <- as.bigq(21,6)
x
# 7 / 2
# Wow ! result is simplified.

y <- as.bigq(5,3)

# addition works !
x + y

# You can even try multiplication, division...
x * y / 13

# and, since May 2012,
x ^ 20
stopifnot(is.bigq(x), is.bigq(x + y),
	  x ^ 20 == as.bigz(7)^20 / 2^20)

# convert to string, double
as.character(x)
as.double(x)

stopifnot( is.na(NA_bigq_) )

# Depict the "S4-class" bigq, i.e., the formal (S4) methods:
if(require("Rmpfr")) # mostly interesting there
  showMethods(class="bigq")

# an  sapply() version that works for big rationals "bigq":
sapplyQ <- function(X, FUN, ...) c_bigq(lapply(X, FUN, ...))

# dummy example showing it works (here):
qq <- as.bigq(1, 1:999)
q1 <- sapplyQ(qq, function(q) q^2)
stopifnot( identical(q1, qq^2) )
}
\keyword{arith}
