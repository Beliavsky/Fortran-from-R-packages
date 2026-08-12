\name{maxsub}
\alias{maxsub}
\alias{maxsub2d}
\title{
Maximal Sum Subarray
}
\description{
Find a subarray with maximal positive sum.
}
\usage{
maxsub(x, inds = TRUE)

maxsub2d(A)
}
\arguments{
  \item{x}{numeric vector.}
  \item{A}{numeric matrix}
  \item{inds}{logical; shall the indices be returned?}
}
\details{
\code{maxsub} finds a contiguous subarray whose sum is maximally positive.
This is sometimes called Kadane's algorithm.
\code{maxsub} will use a very fast version with a running time of
\code{O(n)} where \code{n} is the length of the input vector \code{x}.

\code{maxsub2d} finds a (contiguous) submatrix whose sum of elements is
maximally positive. The approach taken here is to apply the one-dimensional
routine to summed arrays between all rows of \code{A}. This has a run-time
of \code{O(n^3)}, though a run-time of \code{O(n^2 log n)} seems possible
see the reference below.
\code{maxsub2d} can solve a 100-by-100 matrix in a few seconds --
but beware of bigger ones.
}
\value{
Either just a maximal sum, or a list this sum as component \code{sum} plus
the start and end indices as a vector \code{inds}.
}
\references{
Bentley, Jon (1986). ``Programming Pearls'', Column 7.
Addison-Wesley Publ. Co., Reading, MA.

T. Takaoka (2002). Efficient Algorithms for the Maximum Subarray Problem by
Distance Matrix Multiplication. The Australasian Theory Symposion, CATS 2002.
}
\note{
  In special cases, the matrix \code{A} may be sparse or (as in the example
  section) only have one nonzero element in each row and column. Expectation 
  is that there may exists a more efficient (say \code{O(n^2)}) algorithm in 
  these special cases.
}
\author{
  HwB  <hwborchers@googlemail.com>
}
\examples{
##  Find a maximal sum subvector
set.seed(8237)
x <- rnorm(1e6)
system.time(res <- maxsub(x, inds = TRUE))
res

##  Standard example: Find a maximal sum submatrix
A <- matrix(c(0,-2,-7,0, 9,2,-6,2, -4,1,-4,1, -1,8,0,2),
            nrow = 4, ncol = 4, byrow =TRUE)
maxsub2d(A)
# $sum:  15
# $inds: 2 4 1 2 , i.e., rows = 2..4, columns = 1..2

\dontrun{
##  Application to points in the unit square:
set.seed(723)
N <- 50; w <- rnorm(N)
x <- runif(N); y <- runif(N)
clr <- ifelse (w >= 0, "blue", "red")
plot(x, y, pch = 20, col = clr, xlim = c(0, 1), ylim = c(0, 1))

xs <- unique(sort(x)); ns <- length(xs)
X  <- c(0, ((xs[1:(ns-1)] + xs[2:ns])/2), 1)
ys <- unique(sort(y)); ms <- length(ys)
Y  <- c(0, ((ys[1:(ns-1)] + ys[2:ns])/2), 1)
abline(v = X, col = "gray")
abline(h = Y, col = "gray")

A <- matrix(0, N, N)
xi <- findInterval(x, X); yi <- findInterval(y, Y)
for (i in 1:N) A[yi[i], xi[i]] <- w[i]

msr <- maxsub2d(A)
rect(X[msr$inds[3]], Y[msr$inds[1]], X[msr$inds[4]+1], Y[msr$inds[2]+1])
}
}
\keyword{ discrete-optimization }
