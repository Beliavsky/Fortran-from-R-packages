\name{tsp_solver}
%\Rdversion{1.1-4}
\alias{tsp_solver}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Travelling salesperson problem solver}
\description{ Interface to travelling salesperson problem solver. 


Consider an integer linear programming (ILP) formulation of DFJ (Dantzig et al, 1954) used in this research. Let \eqn{G=(V,A)} be a graph with a set \eqn{V} of  \eqn{n} vertices 
and \eqn{A} be a set of arcs or edges. Let \eqn{C=(c_{ij})} be a distance (or cost) matrix associated with \eqn{A}. Elements of the distance matrix \eqn{C,} \eqn{c_{ij},} are positive integers, \eqn{i,j \in V, i \neq j.} 
The TSP focuses on finding a minimum distance circuit (a tour or Hamiltonian circuit) that passes through each vertex once and only once. The DFJ formulation is

 \code{minimize    } \eqn{L =  \sum\limits_{j \neq i}c_{ij}\delta_{ij} ~~~~~~~~~~}   (1)

 \code{subject to    }   \eqn{\sum\limits_{j=1}^n \delta_{ij} = 1, i=1,..., n ~~~~~~~~~~}      (2)

  \code{             } \eqn{~~~~~~~~~~~~~~~~~~~~\sum\limits_{i=1}^n \delta_{ij} = 1, j=1,..., n ~~~~~~~~~~}     (3)

\code{             } \eqn{~~~~~~~~~~~~~~~~~~~~\sum\limits_{i,j\in S}\delta_{ij} \leq |S|-1, S\subset V,  2\leq |S| \leq n-2 ~~~~~~~~~~}    (4)

\code{             } \eqn{~~~~~~~~~~~~~~~~~~~~~ \delta_{ij} \in {0,1}, i,j=1,..., n, i\neq j~~~~~~~~~~~ }      (5)


Constraints (2) and (3) are known as degree constraints indicating that every vertex should be entered and left exactly once correspondingly. Constraints (4) are 
subtour elimination constraints that prevent from forming subtours (several unconnected tours on subsets of less than \eqn{n} vertices), with \eqn{|S|} denoting the number of vertices in \eqn{S}. 

In the DFJ formulation there are \eqn{n(n-1)} unknown binary variables, \eqn{2n} degree constraints and \eqn{2^n -2n-2} subtour elimination constraints. Since the number of subtour elimination constraints increases exponentially, 
solving this problem directly using an integer linear programming code is in general  intractable. However, relaxed versions of the integer linear programming problem where some constraints are initially removed,
and later restored via an iterative process, have been proposed and extensively used.

Here it is proposed to combine heuristics (to get an initial feasible solution) and a linear Diophantine equation (nilde) relaxation to develop a new exact
algorithm that constructs all existing optimal solutions for the TSP in an efficient way.

Below is a brief summary of the proposed algorithm.


Step 1. (Initialization) Solve a corresponding assignment problem to obtain an initial lower bound on the value of the optimal TSP solution. Apply heuristics to obtain an initial upper bound.

Step 2. (Subproblem solution) Given the initial lower bound construct all 0-1 solutions to a linear Diophantine equation introduced by Voinov and Nikulin (1997).

Step 3. (Degree constraints check) Remove solutions  that do not satisfy the degree constraints.  

Step 4. (Subtour elimination) Remove solutions that contain subtours by applying a new simple subtour elimination routine. 
    If there is a solution(s) that contains no subtours, it forms the optimal solution(s): stop. Otherwise, increase the initial lower bound by one and go to step 2. Repeat until the upper bound is reached.
    
  

  
The integer programming formulation of the assignment problem solved in Step 1 of the above algorithm is obtained by relaxing constraints (4), i.e. given by (1) subject to constraints (2), (3) and 
(5).  

For implementing Step 2, solutions of the corresponding subset sum problem should be enumerated. A subset sum problem formulation can be expressed as
\deqn{  a_1s_1+a_2s_2+\ldots +a_ps_p=L, ~~~~~~~~(6)} 
where \eqn{s_i \in \{0,1\},} \eqn{i=1,\ldots,p, p=n(n-1)} is the number of unknown binary variables of the original TSP. \eqn{a_i} are positive integers matching the costs \eqn{c_{ij}} of the cost matrix \eqn{C.}

Voinon and Nikulin (1997)  introduced an algorithm that enumerates all nonnegative integer solutions of equation (6) by using the corresponding generating
function and the binomial theorem. 
All \eqn{0-1} solutions to the equation in (6) can be found by means of the following generating function:
\deqn{
  \Psi_L(z)=\left(z^{a_1}+z^{a_2}+\ldots +z^{a_p} \right)^L= \sum\limits_{k=L\cdot \textrm{min}_i(a_i)}^{k=L\cdot \textrm{max}_i(a_i)}R_k(L,p),
  }
where
\deqn{
   R_k(L,p)= \sum\limits_{s_p=0}^{\textrm{min}\left(1,\left[\frac{L}{a_p}\right]\right)}\sum\limits_{s_{p-1}=0}^{\textrm{min}\left(1,\left[\frac{L-a_ps_p}{a_{p-1}}\right]\right)}
   \ldots \sum\limits_{s_2=0}^{\textrm{min}\left(1,\left[\frac{L-a_ps_p-\ldots -a_3s_3}{a_2}\right]\right)} \frac{L!}{(L-s_1-\ldots -s_p)!s_1!\ldots s_p!}, 
}  
\deqn{ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~(7)
}
\eqn{s_1=\frac{L-a_ps_p-\ldots -a_2s_2}{a_1}} is necessarily either \eqn{0} or \eqn{1}. Otherwise, the equation in (6) does not  have any solutions. The notation \eqn{[x]} denotes the greatest integer part of \eqn{x.}
The right-hand side multiplier in (7) presents the total number of compositions that satisfy the above condition.
If the value of that multiplier is set to \eqn{1}, (7) gives the number of \eqn{0-1} solutions for the equation (6). The solutions, if exist, are written explicitly as
\deqn{
\left\{a_1^{s_1},a_2^{s_2},\ldots, a_p^{s_p} \right\}, ~~~~~~~~(8)
} 
where \eqn{\{s_2,\ldots, s_p\}} are sets of summation indices  in (7), with \eqn{s_1} as specified above. The notation (8) means that in a particular partition 
(a solution of the equation (6)) there are \eqn{s_1} terms equal to \eqn{a_1,} \eqn{s_2} terms of \eqn{a_2} and so on. 



}
\usage{
 tsp_solver(data, labels=NULL,cluster=0,upper_bound=NULL,
          lower_bound=NULL,method="cheapest_insertion",no_go=NULL)

}
\arguments{
   \item{data}{An \code{n x n} matrix of costs/distances of the TSP (with \code{0}'s or \code{NA}s on the main diagonal). Costs/distances of the unconnected edges must be supplied as NA.}
   \item{labels}{An \code{n} vector of optional city labels. If not given, labels are taken from \code{data}.}
   \item{cluster}{Degree constraints can be checked in parallel using parLapply from the parallel package.
      \code{cluster} is either \code{0} (default) for no parallel computing to be used; or \code{1} for one less than the number of cores; or user-supplied cluster on which to do checking. a cluster here can be some cores of a single machine.}
    \item{upper_bound}{A positive integer, an upper bound of the tour length (cost function), if not supplied (default: \code{NULL}) heuristic solution is obtained using 

       \code{TSP::solve_TSP(data,method)}.}
    \item{lower_bound}{A positive integer, a lower possible value of the tour lenght (cost function); 
      if not supplied (default: \code{NULL}), obtained by solving a corresponding assignment problem
      using \code{lpSolve::lp.assign(data)}}
   \item{method}{Heuristic method used in \code{TSP::solve_TSP()} 

    (default: \code{cheapest_insertion})}
   \item{no_go}{A suitably large value used in the distance/cost matrix to make related edges infeasible, 
     if \code{NULL} (default) set to \code{max(data)*10^5}. It can be set to \code{Inf} for \code{TSP()}. However, \code{lpSolve()} is very sensitive to too large values and can result in high values of the \code{lower_bound}.}
}

\value{
  \item{tour}{optimal tour(s).}
  \item{tour_length}{an optimal (minimal) length of the obtained tour(s).}
   \item{coming_solutions}{a list of coming feasible tours obtained within \code{[lower_bound, upper_bound]}.}
  \item{coming_tour_lengths}{a vector of feasible tour length gone within \code{[lower_bound, upper_bound]}.}
  \item{iter}{a number of feasible tour length gone through}
   \item{upper_bound}{an upper bound of the tour length}
   \item{lower_bound}{a lower bound value of the tour lenght}
}

\author{
Vassilly Voinov, Natalya Pya Arnqvist
}
\references{
  Voinov, V. and Nikulin, M. (1997) On a subset sum algorithm and its probabilistic and other applications. In: Advances in combinatorial methods and applications to probability and statistics, Ed. N. Balakrishnan, \enc{Birkhäuser}{Birkhauser}, Boston, 153-163.

 Dantzig, G., Fulkerson, R. and Johnson, S. (1954) Solution of a large-scale traveling-salesman problem.
Journal of the operations research society of America , 2(4):393-410.

}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
    \code{\link{nilde-package}}, \code{\link{get.partitions}}, \code{\link{get.knapsack}}, \code{\link{get.subsetsum}}     
}

\examples{
\dontrun{
## some examples...
library(nilde)
set.seed(3)
d <- matrix(sample(1:100,25,replace=TRUE),5,5)
diag(d) <-NA # although no need to specify as the code assumes NAs by default
g <- tsp_solver(d)
g
}

}
\keyword{optimize}% at least one, from doc/KEYWORDS



