\name{maxempty}
\alias{maxempty}
\title{
  Maximally Empty Rectangle Problem
}
\description{
  Find the largest/maximal empty rectangle, i.e. with largest area, not 
  containing given points.
}
\usage{
maxempty(x, y, ax = c(0, 1), ay = c(0, 1))
}
\arguments{
  \item{x, y}{coordinates of points to be avoided.}
  \item{ax, ay}{left and right resp. lower and upper constraints.}
}
\details{
  Find the largest or maximal empty two-dimensional rectangle in a rectangular
  area. The edges of this rectangle have to be parallel to the edges of the
  enclosing rectangle (and parallel to the coordinate axes). `Empty' means
  that none of the points given are contained in the interior of the found
  rectangle.
}
\value{
  List with \code{area} and \code{rect} the rectangle as a vector usable for
  the \code{rect} graphics function.
}
\note{
  The algorithm has a run-time of \code{O(n^2)} while there are run-times of
  \code{O(n*log(n))} reported in the literature, utilizing a more complex
  data structure. I don't know of any comparable algorithms for the largest 
  empty circle problem.
}
\author{
  HwB  email: <hwborchers@googlemail.com>
}
\references{
  B. Chazelle, R. L. Drysdale, and D. T. Lee (1986). Computing the Largest 
  Empty Rectangle. SIAM Journal of Computing, Vol. 15(1), pp. 300--315.

  A. Naamad, D. T. Lee, and W.-L. Hsu (1984). On the Maximum Empty Rectangle 
  Problem. Discrete Applied Mathematics, Vol. 8, pp. 267--277.
}
\seealso{
  \code{Hmisc::largest.empty} with a Fortran implementation of this code.
}
\examples{
N <- 100; set.seed(8237)
x <- runif(N); y <- runif(N)
R <- maxempty(x, y, c(0,1), c(0,1))
R
# $area
# [1] 0.08238793
# $rect
# [1] 0.7023670 0.1797339 0.8175771 0.8948442

\dontrun{
plot(x, y, pch="+", xlim=c(0,1), ylim=c(0,1), col="darkgray",
      main = "Maximally empty rectangle")
rect(0, 0, 1, 1, border = "red", lwd = 1, lty = "dashed")
do.call(rect, as.list(R$rect))
grid()}
}
\keyword{ discrete-optimization }
