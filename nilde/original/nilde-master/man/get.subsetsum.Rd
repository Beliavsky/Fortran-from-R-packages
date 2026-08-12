\name{get.subsetsum}
%\Rdversion{1.0}
\alias{get.subsetsum}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Enumeration of all existing 0-1 and bounded solutions of a subset sum problem}
\description{By default this function solves the following 0-1 subset sum problem.
Given the set of positive integers \code{(a_1, a_2, ..., a_l)} and a positive integer \code{n}, find all non-empty subsets that sum to \code{n}, so that each of the integers \code{a_i} either appears in the subset or it does not, and the total number of summands should not exceed \code{M}, \code{M <= n}.

The bounded subset sum problem has restrictions on the number of times (bounds) \code{a_i} can turn up in the subset.
  
 The algorithm is based on a generating function of Hardy and Littlewood used by Voinov and Nikulin (1997). 
}
\usage{
 get.subsetsum(a,n,M=NULL,problem="subsetsum01",bounds=NULL)

}
\arguments{
   \item{a}{An \code{l}-vector of positive integers with \code{l>= 2}.}
   \item{n}{A positive integer.}
    \item{M}{A positive integer, the maximum number of summands, \code{M <= n}} 
    \item{problem}{one of the two problems to be solved: "subsetsum01" (default) for a 0-1 subset sum problem, 
         or "bsubsetsum" a bounded subset sum problem.}
    \item{bounds}{An \code{l}-vector of positive integers, bounds for \eqn{s_i}, i.e. \eqn{0 <= s_i <= b_i}}
}

\value{
  \item{p.n}{total number of solutions obtained.}
  \item{solutions}{a matrix with each column presenting a solution for \code{n}.}
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
    \code{\link{nilde-package}}, \code{\link{get.partitions}}, \code{\link{get.knapsack}}, \code{\link{nlde}}     
}

\examples{
## some examples...
b1 <- get.subsetsum(a=c(41,34,21,20,8,7,7,4,3,3),M=10,n=50,problem="subsetsum01")
b1
colSums(b1$solutions*c(41,34,21,20,8,7,7,4,3,3))

b2 <- get.subsetsum(a=c(111:120),M=10,n=485,problem="subsetsum01") ## no solutions
b2

b3 <- get.subsetsum(a=c(30,29,32,31,33),M=5,n=91,problem="subsetsum01")
b3
colSums(b3$solutions*c(30,29,32,31,33))
get.subsetsum(a=c(30,29,32,31,33),M=5,n=91,problem="bsubsetsum",bounds=c(1,1,1,1,1))

b4 <- get.subsetsum(a=c(30,29,32,31,33),M=5,n=91,problem="bsubsetsum",
                 bounds=c(1,2,1,3,4))
b4
colSums(b4$solutions*c(30,29,32,31,33))
}
\keyword{optimize}% at least one, from doc/KEYWORDS



