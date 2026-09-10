\name{asNumeric}
\alias{asNumeric}
\alias{asNumeric-methods}
\alias{asNumeric,ANY-method}
\alias{asNumeric,bigq-method}
\alias{asNumeric,bigz-method}
\title{Coerce to 'numeric', not Loosing Dimensions}
\description{
  a number-like object is coerced to type (\link{typeof})
  \code{"numeric"}, keeping \code{\link{dim}} (and maybe
  \code{\link{dimnames}}) when present.
}
\usage{
asNumeric(x)
}
\arguments{
  \item{x}{a \dQuote{number-like} object, e.g., big integer
    (\code{\link{bigz}}), or \code{\link[Rmpfr]{mpfr}}, notably
    including matrices and arrays of such numbers.}
}
\section{Methods}{
  \describe{
    \item{\code{signature(x = "ANY")}}{the default method, which is the
      identity for \code{\link{numeric}} array.}
    \item{\code{signature(x = "bigq")}}{the method for big rationals.}
    \item{\code{signature(x = "bigq")}}{the method for big integers.}
  }
  Note that package \pkg{Rmpfr} provides methods for its own number-like objects.
  % -->  ~/R/Pkgs/Rmpfr/man/asNumeric-methods.Rd
}
%% \details{
%% %%  ~~ If necessary, more details than the description above ~~
%% }
\value{
  an \R object of type (\code{\link{typeof}}) \code{"numeric"}, a \code{\link{matrix}}
  or \code{\link{array}} if \code{x} had non-NULL dimension \code{\link{dim}()}.
}
\author{Martin Maechler}

\seealso{
 \code{\link{as.numeric}} coerces to both \code{"numeric"} and to a
 \code{\link{vector}}, whereas \code{asNumeric()} should keep
 \code{\link{dim}} (and other) attributes.
}
\examples{
m <- matrix(1:6, 2,3)
stopifnot(identical(m, asNumeric(m)))# remains matrix

(M <- as.bigz(m) / 5) ##-> "bigq" matrix
asNumeric(M) # numeric matrix
stopifnot(all.equal(asNumeric(M), m/5))
}
\keyword{arith}
\keyword{methods}

