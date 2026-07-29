\name{returns}
\alias{returns}
\alias{returns_qrmtools}
\title{Computing Returns and Inverse Transformation}
\description{
  Compute log-returns, simple returns and basic differences (or the
  inverse operations) from given data.
}
\usage{
returns(x, method = c("logarithmic", "simple", "diff"), inverse = FALSE,
        start, start.date)
returns_qrmtools(x, method = c("logarithmic", "simple", "diff"),
                 inverse = FALSE, start, start.date)
}
\arguments{
  \item{x}{matrix or vector (possibly a \code{\link[xts]{xts}} object) to be turned
    into returns (if \code{inverse = FALSE})
    or returns to be turned into the original data (if \code{inverse =
      TRUE}).}
  \item{method}{\code{\link{character}} string indicating the method to
    be used (log-returns (logarithmic changes), simple returns (relative
    changes), or basic differences). Note that this can also be a vector
    of such methods of length equal to the number of columns of \code{x}.}
  \item{inverse}{\code{\link{logical}} indicating whether the inverse
    transformation (data from given returns) shall be computed (if
    \code{TRUE}, this requires \code{start} to be specified).}
  \item{start}{if \code{inverse = TRUE}, the last available value
    of the time series to be constructed from the given returns \code{x}.}
  \item{start.date}{\code{\link{character}} or \code{\link{Date}} object
    to be used as the date corresponding to the value \code{start};
    currently only used for \code{\link[xts]{xts}} objects.}
}
\value{
  \code{\link{vector}} or \code{\link{matrix}} with the same number of
  columns as \code{x} just one row less if \code{inverse = FALSE}
  or one row more if \code{inverse = TRUE}.
}
\details{
  If \code{inverse = FALSE} and \code{x} is an \code{\link[xts]{xts}}
  object, the returned object is an \code{\link[xts]{xts}}, too.

  Note that the \R package \pkg{timeSeries} also contains a function
  \code{returns()} (and hence the order in which \pkg{timeSeries} and
  \pkg{qrmtools} are loaded matters to get the right \code{returns()}).
  For this reason, \code{returns_qrmtools()} is an alias for
  \code{returns()} from \pkg{qrmtools}.
}
\author{Marius Hofert}
\examples{
## Generate two paths of a geometric Brownian motion
S0 <- 10 # current stock price S_0
r <- 0.01 # risk-free annual interest rate
sig <- 0.2 # (constant) annual volatility
T <- 2 # maturity in years
N <- 250 # business days per year
t <- 1:(N*T) # time points to be sampled
npath <- 2 # number of paths
set.seed(271) # for reproducibility
S <- replicate(npath, S0 * exp(cumsum(rnorm(N*T, # sample paths of S_t
                                            mean = (r-sig^2/2)/N,
                                            sd = sqrt((sig^2)/N))))) # (N*T, npath)

## Turn into xts objects
library(xts)
sdate <- as.Date("2000-05-02") # start date
S.  <- as.xts(S, order.by = seq(sdate, length.out = N*T, by = "1 week"))
plot(S.[,1], main = "Stock 1")
plot(S.[,2], main = "Stock 2")


### Log-returns ################################################################

## Based on S[,1]
X <- returns(S[,1]) # build log-returns (one element less than S)
Y <- returns(X, inverse = TRUE, start = S[1,1]) # transform back
stopifnot(all.equal(Y, S[,1]))

## Based on S
X <- returns(S) # build log-returns (one element less than S)
Y <- returns(X, inverse = TRUE, start = S[1,]) # transform back
stopifnot(all.equal(Y, S))

## Based on S.[,1]
X <- returns(S.[,1])
Y <- returns(X, inverse = TRUE, start = S.[1,1], start.date = sdate)
stopifnot(all.equal(Y, S.[,1], check.attributes = FALSE))

## Based on S.
X <- returns(S.)
Y <- returns(X, inverse = TRUE, start = S.[1], start.date = sdate)
stopifnot(all.equal(Y, S., check.attributes = FALSE))

## Sign-adjusted (negative) log-returns
X <- -returns(S) # build -log-returns
Y <- returns(-X, inverse = TRUE, start = S[1,]) # transform back
stopifnot(all.equal(Y, S))


### Simple returns #############################################################

## Simple returns based on S
X <- returns(S, method = "simple")
Y <- returns(X, method = "simple", inverse = TRUE, start = S[1,])
stopifnot(all.equal(Y, S))

## Simple returns based on S.
X <- returns(S., method = "simple")
Y <- returns(X, method = "simple", inverse = TRUE, start = S.[1,],
             start.date = sdate)
stopifnot(all.equal(Y, S., check.attributes = FALSE))

## Sign-adjusted (negative) simple returns
X <- -returns(S, method = "simple")
Y <- returns(-X, method = "simple", inverse = TRUE, start = S[1,])
stopifnot(all.equal(Y, S))


### Basic differences ##########################################################

## Basic differences based on S
X <- returns(S, method = "diff")
Y <- returns(X, method = "diff", inverse = TRUE, start = S[1,])
stopifnot(all.equal(Y, S))

## Basic differences based on S.
X <- returns(S., method = "diff")
Y <- returns(X, method = "diff", inverse = TRUE, start = S.[1,],
             start.date = sdate)
stopifnot(all.equal(Y, S., check.attributes = FALSE))

## Sign-adjusted (negative) basic differences
X <- -returns(S, method = "diff")
Y <- returns(-X, method = "diff", inverse = TRUE, start = S[1,])
stopifnot(all.equal(Y, S))


### Vector-case of 'method' ####################################################

X <- returns(S., method = c("logarithmic", "diff"))
Y <- returns(X, method = c("logarithmic", "diff"), inverse = TRUE, start = S.[1,],
             start.date = sdate)
stopifnot(all.equal(Y, S., check.attributes = FALSE))
}
\keyword{utilities}