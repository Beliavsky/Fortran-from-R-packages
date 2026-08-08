\name{next_subset}
\alias{next_subset}
\title{
  Select Next Subset
}
\description{
Select next subset of size k from a set of size n.
}
\usage{
next_subset(a, n, k)
}
\arguments{
  \item{a}{
    a numeric vector (integers)
}
  \item{n}{
    an integer: the size of the set to choose from
}
  \item{k}{
    an integer: the subset size
}
}
\details{

  Given a subset \emph{a} of size \emph{k} taken from a set
  of size \emph{n}, compute the next subset by alphabetical
  order.

  Uses algorithm NEXKSB of Nijenhuis and Wilf (1975).
}
\value{
   a numeric vector (the next subset) or \code{\link{NULL}}
   (when there is no next subset)
}
\references{
  Nijenhuis, A. and Wilf, H. S. (1975)
  \emph{Combinatorial Algorithms for Computers and Calculators}.
  Academic Press.
}
\author{
  Enrico Schumann
}
\seealso{

  \code{\link{choose}} computes the number of combinations\cr
  \code{\link{combn}} creates all combinations\cr
  \code{\link{expand.grid}}
}
\examples{
n <- 4
k <- 2
t(combn(n, k))
##      [,1] [,2]
## [1,]    1    2
## [2,]    1    3
## [3,]    1    4
## [4,]    2    3
## [5,]    2    4
## [6,]    3    4

a <- 1:k
print(a)
while (!is.null(a))
    print(a <- next_subset(a, n = n, k = k))
## [1] 1 2
## [1] 1 3
## [1] 1 4
## [1] 2 3
## [1] 2 4
## [1] 3 4
}
