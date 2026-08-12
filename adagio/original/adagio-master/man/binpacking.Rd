\name{bpp_approx}
\alias{bpp_approx}
\title{
  Approximate Bin Packing
}
\description{
  Solves the Bin Packing problem approximately.
}
\usage{
bpp_approx(S, cap, method = c("firstfit", "bestfit", "worstfit"))
}
\arguments{
  \item{S}{vector of weights (or sizes) of items.}
  \item{cap}{same capacity for all the bins.}
  \item{method}{which approximate method to use.}
}
\details{
  Solves approximately the Bin Packing problem for numeric weights and bins,
  all having the same volume.

  Possible methods are "firstfit", "bestfit", and "worstfit". "firstfit" tries
  to place each item as early as possible, "bestfit" such that the remaining
  space in the bin is as small as possible, and "worstfit" such that the
  remaining space is as big as possible.

  Best results are achieved with the "bestfit" method. "firstfit" may be a
  reasonable alternative. For smaller and medium-sized data the approximate
  results will come quite close to the exact solution, see the examples.

  In general, the results are much better if the items in \code{S} are sorted
  decreasingly. If they are not, an immediate warning is issued.
}
\value{
  A list of the following components:
  \item{nbins}{minimum number of bins.}
  \item{xbins}{index of the bin each item is assigned to.}
  \item{sbins}{sum of item sizes in each bin.}
  \item{filled}{total volume filled in the bins (as percentage).}
}
\references{
  Silvano Martello. "Bin packing problems". In: 23rd Belgian Mathematical 
  Optimization Workshop, La-Roche-en-Ardennes 2019.
}
\author{
  Hans W. Borchers
}
\note{
  The Bin Packing problem can be solved as a Linear Program. The formulation
  is a bit tricky, and it turned out 'lpSolve' does not solve medium-sized
  problems in acceptable time. (Tests with 'Rglpk' will follow.)
}
\seealso{
  Function \code{binpacking} in package 'knapsack' (on R-Forge).
}
\examples{
## (1)
S <- c(50, 3, 48, 53, 53, 4, 3, 41, 23, 20, 52, 49)
cap <- 100
bpp_approx(S, cap, method = "bestfit")
## exact    -- $nbins 4, filled 99.75 %
## firstfit -- $nbins 6, filled 66.5  %
## bestfit  -- $nbins 5, filled 79.8  %
## ! when decreasingly sorted, 'bestfit' with nbins = 4

## (2)
S <- c(100,99,89,88,87,75,67,65,65,57,57,49,47,31,27,18,13,9,8,1)
cap <- 100
bpp_approx(S, cap, method = "firstfit")
# firstfit: 12 bins; exact: 12 bins

\dontrun{
## (3)
S <-  c(99,99,96,96,92,92,91,88,87,86,
        85,76,74,72,69,67,67,62,61,56,
        52,51,49,46,44,42,40,40,33,33,
        30,30,29,28,28,27,25,24,23,22,
        21,20,17,14,13,11,10, 7, 7, 3)
cap <- 100
bpp_approx(S, cap)
# exact: 25; firstfit: 25; bestfit: 25 nbins

## (4)
# 20 no.s in 1..100, capacity 100
set.seed(7013)
S <- sample(1:100, 20, replace = TRUE)
cap <- 100
bpp_approx(sort(S, decreasing = TRUE), cap, method = "bestfit")
# exact: 12 bins; firstfit and bestfit: 13; worstfit: 14 bins
}
}
\keyword{ discrete-optimization }
