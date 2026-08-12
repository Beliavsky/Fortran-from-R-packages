\name{get.knapsack}
%\Rdversion{1.0}
\alias{get.knapsack}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Enumeration of all existing nonnegative integer  solutions for unbounded, bounded and 0-1 knapsack and subset sum  problems}
\description{This function solves the unbounded, bounded and 0-1 knapsack problems.

The unbounded knapsack problem can be written as follows.

 \code{maximize   } \eqn{c_1s_1 +c_2s_2 +...+ c_ls_l }

 \code{subject to } \eqn{a_1s_1 +a_2s_2 +...+ a_ls_l <= n,}

 \code{             }\eqn{s_i >= 0, integers.}

The bounded knapsack problem has additional constraints, \eqn{0 <= s_i <= b_i, i=1,...,l, b_i <= [n/a_i].} The 0-1 knapsack problem arises when 
\eqn{s_i= 0} or \eqn{1, i=1,...,l}.  

 The algorithm is based on a generating function of Hardy and Littlewood used by Voinov and Nikulin (1997). Subset sum problems are particular
cases of knapsack problems when vectors of \code{weights}, \eqn{(a_1,...,a_l)}, and \code{objectives}, \eqn{(c_1,...,c_l)}, are equal.
}
\usage{
 get.knapsack(objective,a,n,problem="uknap",bounds=NULL)

}
\arguments{
   \item{objective}{A vector of coefficients (values of each item in the knapsack) of the objective function to be maximized when solving knapsack problem.}
   \item{a}{An \code{l}-vector of weights of each item in a knapsack, with \code{l>= 2}.}
   \item{n}{a maximal possible capacity of the knapsack.} 
    \item{problem}{one of the following names of the problems to be solved: 
         "uknap" (default) for an unbounded knapsack problem,
          "knap01" for a 0-1 knapsack problem, and
         "bknap" for a bounded knapsack problem.}
    \item{bounds}{An l-vector of positive integers, bounds of \eqn{s_i}, i.e. \eqn{0 <= s_i <= b_i}.}
}

\value{
  \item{p.n}{total number of solutions obtained.}
  \item{solutions}{a matrix with each column representing a solution of \code{n}.}
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
    \code{\link{nilde-package}}, \code{\link{get.partitions}}, \code{\link{get.subsetsum}}, \code{\link{nlde}}    
}

\examples{
## some examples...
b1 <- get.knapsack(objective=c(200:206),a=c(100:106),n=999,problem="uknap")
b1

b2 <- get.knapsack(objective=c(41,34,21,20,8,7,7,4,3,3),a=c(41,34,21,20,8,7,7,4,3,3),
         n=50, problem="bknap", bounds=rep(2,10))
b2
colSums(b2$solutions*c(41,34,21,20,8,7,7,4,3,3))

b3 <- get.knapsack(objective=c(4,3,3),a=c(3,2,2),n=4,problem="bknap",bounds=c(2,2,2))
b3
## get maximum value of the objective function...
colSums(b3$solutions*c(4,3,3)) 
## checking constraint...
colSums(b3$solutions*c(3,2,2))   

b4 <- get.knapsack(objective=c(4,3,3),a=c(3,2,2),n=4,problem="knap01")
b4
## get maximum value of the objective function...
colSums(b4$solutions*c(4,3,3)) 
## checking constraint...
colSums(b4$solutions*c(3,2,2))

\dontrun{
b5 <- get.knapsack(a=c(100:106),n=2999,objective=c(200:206),problem="uknap")
b5$p.n ## total number of solutions
options(max.print=5E5)
print(b5)

  }
}
\keyword{optimize}% at least one, from doc/KEYWORDS




