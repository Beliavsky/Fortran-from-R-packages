\name{mpfr}
\title{Exported function for mpfr use}
\alias{.as.bigz}
\alias{..as.bigz}
\alias{.as.char.bigz}
\alias{.sub.bigq}

\description{
	Theses hidden function are provided for mpfr use. Use theses function with care.

}
\usage{
.as.bigz(a, mod=NA)
}
\arguments{
 \item{a}{either \code{\link{integer}}, \code{\link{numeric}} (i.e.,
    \code{\link{double}}) or \code{\link{character}} vector.

    If character: the strings either start with \code{0x} for hexadecimal, \code{0b} for
    binary, \code{0} for octal, or without a \code{0*} prefix for decimal values.
    Formatting errors are signalled as with \code{\link{stop}}.}
  \item{mod}{an integer, numeric, string or bigz of the internal modulus, see below.}

 }
\references{The GNU MP Library, see \url{https://gmplib.org}}
% \details{
%   ~~ If necessary, more details than the description above ~~
% }
\value{
  An \R object of (S3) class \code{"bigz"}, representing the argument
  (\code{x} or \code{a}).
}
\examples{
.as.bigz(1)
}
\keyword{arith}
