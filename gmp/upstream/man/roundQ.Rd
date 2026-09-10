\name{roundQ}
\title{Rounding Big Rationals ("bigq") to Decimals}
\alias{roundQ}
\alias{round0}
\alias{round.bigq}
\description{
  Rounding big rationals (of class \code{"bigq"}, see \code{\link{as.bigq}()})
  to decimal \code{digits} is strictly based on a (optionally choosable)
  definition of rounding to integer, i.e., \code{digits = 0}, the default
  method of which we provide as \code{round0()}.

  The users typically just call \code{round(x, digits)} as elsewhere, and
  the \code{round()} method will call \code{round(x, digits, round0=round0)}.
}
\usage{
round0(x)

roundQ(x, digits = 0, r0 = round0)

\S3method{round}{bigq}(x, digits = 0)
}
\arguments{
  \item{x}{vector of big rationals, i.e., of \code{\link{class}} \code{"bigq"}.}
  \item{digits}{integer number of decimal digits to round to.}
  \item{r0}{a \code{\link{function}} of one argument which implements a
    version of \code{\link{round}(x, digits=0)}.   The default for
    \code{roundQ()} is to use our \code{round0()} which implements
    \dQuote{round to even}, as base \R's \code{\link{round}}.}
}
\value{
  \code{round0()} returns a vector of big integers, i.e., \code{"bigz"} classed.

  \code{roundQ(x, digits, round0)} returns a vector of big rationals,
  \code{"bigq"}, as \code{x}.

  \code{round.bigq} is very simply defined as
  \code{function(x, digits)  roundQ(x, digits)} .
}
\references{
  The vignette \dQuote{\emph{Exact Decimal Rounding via Rationals}} from
  CRAN package \CRANpkg{round},
  % not yet
  % \url{https://CRAN.R-project.org/package=round/vignettes/rationalRound.html}

  Wikipedia, Rounding, notably "Round half to even":
  \url{https://en.wikipedia.org/wiki/Rounding#Round_half_to_even}
}
\author{Martin Maechler, ETH Zurich}
\seealso{
  \code{\link{round}} for (double precision) numbers in base \R;
  \code{\link[round]{roundX}} from CRAN package \CRANpkg{round}.
}
\examples{
qq <- as.bigq((-21:31), 10)
noquote(cbind(as.character(qq), asNumeric(qq)))
round0(qq) # Big Integer ("bigz")
## corresponds to R's own "round to even" :
stopifnot(round0(qq) == round(asNumeric(qq)))
round(qq) # == round(qq, 0): the same as round0(qq) *but* Big Rational ("bigq")

halfs <- as.bigq(1,2) + -5:12
\dontshow{q <- c(halfs, qq)
stopifnot(round0(q) == round(q)) ;  rm(q)
}

\dontshow{if(FALSE)}% do not create it in user's globalenv
## round0() is simply
round0 <- function (x) {
    nU <- as.bigz.bigq(xU <- x + as.bigq(1, 2)) # traditional round: .5 rounded up
    if(any(I <- is.whole.bigq(xU))) { # I <==>  x == <n>.5 : "hard case"
        I[I] <- .mod.bigz(nU[I], 2L) == 1L # rounded up is odd  ==> round *down*
        nU[I] <- nU[I] - 1L
    }
    nU
}

## 's' for simple: rounding as you learned in school:
round0s <- function(x) as.bigz.bigq(x + as.bigq(1, 2))

cbind(halfs, round0s(halfs), round0(halfs))

\dontshow{if(FALSE)}
## roundQ() is simply
roundQ <- function(x, digits = 0, r0 = round0) {
    ## round(x * 10^d) / 10^d --  vectorizing in both (x, digits)
    p10 <- as.bigz(10) ^ digits # class: if(all(digits >= 0)) "bigz" else "bigq"
    r0(x * p10) / p10
}
}
\keyword{arith}
\concept{Rounding}
