\name{binomQ}
\title{Exact Rational Binomial Probabilities}
\alias{dbinomQ}
%NOT YET \alias{pbinomQ}
\description{
  Compute \emph{exact} binomial probabilities using (big integer and) big
  rational arithmetic.
}
\usage{
dbinomQ(x, size, prob, log = FALSE)
}
%% pbinomQ(q, size, prob, lower.tail = TRUE, log.p = FALSE)
%
\arguments{
  \item{x, size}{integer or big integer (\code{"\link{bigz}"}), will be
    passed to \code{\link{chooseZ}()}.}
  \item{prob}{the probability; should be big rational
    (\code{"\link{bigq}"}); if not it is coerced with a warning.}
  \item{log}{logical; must be \code{FALSE} on purpose.  Use
    \code{log(Rmpfr::\link[Rmpfr]{mpfr}(dbinomQ(..), precB))} for the
    logarithm of such big rational numbers.}
}
\value{a big rational (\code{"\link{bigq}"}) of the \code{\link{length}} of
  (recycled) \code{x+size+prob}.
}
\author{Martin Maechler}
\seealso{
  \code{\link{chooseZ}}; \R's (\pkg{stats} package) \code{\link{dbinom}()}.
}
\examples{
dbinomQ(0:8,8, as.bigq(1,2))
##  1/256  1/32   7/64   7/32   35/128 7/32   7/64   1/32   1/256

ph16. <- dbinomQ(0:16, size=16, prob = 1/2)  # innocous warning
ph16  <- dbinomQ(0:16, size=16, prob = as.bigq(1,2))
ph16.75 <- dbinomQ(0:16, size=16, prob = as.bigq(3,4))
ph8.75  <- dbinomQ(0:8, 8, as.bigq(3,4))
stopifnot(exprs = {
   dbinomQ(0:8,8, as.bigq(1,2)) * 2^8 == choose(8, 0:8)
   identical(ph8.75, chooseZ(8,0:8) * 3^(0:8) / 4^8)
   all.equal(ph8.75, choose (8,0:8) * 3^(0:8) / 4^8, tol=1e-15) # see exactly equal
   identical(ph16, ph16.)
   identical(ph16,
            dbinomQ(0:16, size=16, prob = as.bigz(1)/2))
   all.equal(dbinom(0:16, 16, prob=1/2), asNumeric(ph16),    tol=1e-15)
   all.equal(dbinom(0:16, 16, prob=3/4), asNumeric(ph16.75), tol=1e-15)
})
}
