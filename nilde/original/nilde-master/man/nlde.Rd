\name{nlde}
%\Rdversion{1.0}
\alias{nlde}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Enumeration of all existing nonnegative integer solutions of a linear Diophantine equation}
\description{This function enumerates nonnegative integer solutions of a linear Diophantine equation (NLDE):
  \deqn{a_1s_1 +a_2s_2 +...+ a_ls_l =n,}
where \eqn{a_1 <= a_2 <= ... <= a_l,} \eqn{a_i > 0,} \eqn{n > 0,} \eqn{s_i >= 0,} \eqn{i=1,2,...,l,} and all variables involved are integers.

 The algorithm is based on a generating function of Hardy and Littlewood used by Voinov and Nikulin (1997). 
}
\usage{
 nlde(a, n, M=NULL, at.most=TRUE, option=0)

}
\arguments{
   \item{a}{An \code{l}-vector of positive integers (coefficients of the left-hand-side of NLDE) with \code{l>= 2}.}
   \item{n}{A positive integer which is to be partitioned.}
   \item{M}{A positive integer, the number of parts of \code{n}, \code{M <= n}.}
   \item{at.most}{If \code{TRUE} partitioning of \code{n} into at most \code{M} parts, if \code{FALSE} partitioning on exactly \code{M} parts.} 
   \item{option}{When set to \code{1} (or any positive number) finds only 0-1 solutions of the linear Diophantine equation.
                 When set to \code{2} (or any positive number > 1) finds 0-1 solutions of the linear Diophantine inequality.}
}
\value{
  \item{p.n}{total number of partitions obtained.}
  \item{solutions}{a matrix with each column forming a partition of \code{n}.}
}

\author{
Vassilly Voinov, Natalya Pya Arnqvist, Yevgeniy Voinov
}
\references{
  Voinov, V. and Nikulin, M. (1997) On a subset sum algorithm and its probabilistic and other applications. In: Advances in combinatorial methods and applications to probability and statistics, Ed. N. Balakrishnan, \enc{Birkhäuser}{Birkhauser}, Boston, 153-163.

  Hardy, G.H. and Littlewood, J.E. (1966) Collected Papers of G.H. Hardy, Including Joint Papers with J.E. Littlewood and Others. Clarendon Press, Oxford.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{nilde-package}}, \code{\link{get.partitions}}, \code{\link{get.subsetsum}}, \code{\link{get.knapsack}}     
}

\examples{
## some examples...
## example 1...
nlde(a=c(3,2,5,16),n=18,at.most=TRUE)
b1 <- nlde(a=c(3,2,5,16),n=18,M=6,at.most=FALSE)
b1
## checking M, the number of parts that n=18 has been partitioned into...
colSums(b1$solutions)
## checking the value of n...
colSums(b1$solutions*c(3,2,5,16))

## example 2:  solving 0-1 nlde ...
b2 <- nlde(a=c(3,2,5,16),n=18,M=6,option=1)
b2
colSums(b2$solutions*c(3,2,5,16))

## example 3...
b3 <- nlde(c(15,21),261)
b3 
## checking M, the number of parts that n has been partitioned into...
colSums(b3$solutions)
## checking the value of n...
colSums(b3$solutions*c(15,21))

## example 4... 
nlde(c(5,6),19) ## no solutions

## example 5: solving 0-1 inequality...
b4 <- nlde(a=c(70,60,50,33,33,33,11,7,3),n=100,at.most=TRUE,option=2)


}
\keyword{optimize}% at least one, from doc/KEYWORDS

