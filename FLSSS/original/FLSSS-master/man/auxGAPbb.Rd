\name{auxGAPbb}
\alias{auxGAPbb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multithreaded generalized assignment problem solver via branch and bound
}
\description{
Multithreaded exact solver for the generalized assignment problem via decomposition to binary knapsack problems (branch), and Lagrangian relaxation (bound).
}
\usage{
auxGAPbb(
  cost,
  profitOrLoss,
  budget,
  maxCore = 7,
  tlimit = 60,
  ub = "MT",
  greedyBranching = TRUE,
  optim = "max",
  multhreadOn = "nodes",
  threadLoad = 32
)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cost}{
A numeric matrix. Dimensionality = N(agents) \code{x} N(tasks).
}
  \item{profitOrLoss}{
A numeric matrix of the same dimensionality of \code{cost}. Profit for maximum GAP. Loss for minimum GAP.
}
  \item{budget}{
A numeric vector. Size = N(agents).
}
  \item{maxCore}{
Maximal threads to invoke. No greater than the number of logical CPUs on machine.
}
  \item{tlimit}{
Return the best exsisting solution in \code{tlimit} seconds.
}
  \item{ub}{
Upper bound function. "MT" or "HS". See \code{auxKnapsack01bb()}.
}
  \item{greedyBranching}{
If \code{TRUE}, branch and bound in a greedy manner. See Details.
}
  \item{optim}{
A string. \code{optim = "max"} (\code{"min"}) solves the maximum (minimum) GAP.
}
\item{multhreadOn}{
A string. The default \code{multhreadOn = "nodes"} multithreads over branching nodes. Internally, a single-threaded miner runs at first and stops once there are no less than \code{maxCore} \code{*} \code{threadLoad} latent trees stored in stack. The miner realizes those branches and distribute them to threads. Threads work on the subproblems sequentially and update the optimum supervised by a light mutex lock.

Other values of \code{multhreadOn} assign threads to knapsack problems at each branching node, which is a historical remain and should always be avoided. It has overwhelming overheads because knapsack problems at those nodes are typically trivial.
}
\item{threadLoad}{
An integer. Each thread is loaded with \code{threadLoad} sub-problems on average.
}
}
\details{
A popular library of GAP instances can be found here: https://github.com/WhateverLiu/gapInstances.

This algorithm is based on a foundational paper by Ross and Soland (1975) and is carefully engineered towards speed. Implementation highlights include (i) decomposition for multithreading; (ii) a new branching method (\code{greedyBranching}) that pushes all candidate branching variables at each node into stack instead of pushing only those that have the highest desirabilities and would not affect the subsequent branching after being pushed; (iii) the return of current best solutions in time; (iv) the capability of taking real costs and profits. \code{greedyBranching} may considerably lower the number of nodes having the same series of knapsack problems to solve, thus accelerate the convergence speed.
}
\value{
A list of 5:

\item{totalProfitOrLoss }{Total profit or loss generated from the assignment.}

\item{agentCost }{A numeric vector of total costs for each agent.}

\item{assignment }{An integer vector. \code{assignment[i]} indexes the agent assigned to the \code{i}th task.}

\item{nodes }{The number of branching nodes generated in mining.}

\item{bkpSolved }{The number of binary knapsack problems solved in mining.}
}

\note{
The C++ implementation is fully independent and borrows no code from any commercial or open source.
}

\examples{
# =============================================================================
# Data source: http://people.brunel.ac.uk/~mastjjb/jeb/orlib/gapinfo.html,
# gap1 c515-1, 5 agents 15 tasks. Parsed instances from the library can be
# found here: https://github.com/WhateverLiu/gapInstances
# =============================================================================
profit = c(17,21,22,18,24,15,20,18,19,18,16,22,24,24,16,23,16,21,16,17,16,19,
           25,18,21,17,15,25,17,24,16,20,16,25,24,16,17,19,19,18,20,16,17,21,
           24,19,19,22,22,20,16,19,17,21,19,25,23,25,25,25,18,19,15,15,21,25,
           16,16,23,15,22,17,19,22,24)
profit = t(matrix(profit, ncol = 5))
cost = c(8,15,14,23,8,16,8,25,9,17,25,15,10,8,24,15,7,23,22,11,11,12,10,17,16,
         7,16,10,18,22,21,20,6,22,24,10,24,9,21,14,11,14,11,19,16,20,11,8,14,
         9,5,6,19,19,7,6,6,13,9,18,8,13,13,13,10,20,25,16,16,17,10,10,5,12,23)
cost = t(matrix(cost, ncol = 5))
budget = c(36, 34, 38, 27, 33)


sol = FLSSS::auxGAPbb(cost, profit, budget, maxCore = 2, tlimit = 4,
                       ub = "MT", greedyBranching = TRUE, optim = "max")




# =============================================================================
# Data source: http://support.sas.com/documentation/cdl/en/ormpug/65554/HTML
# /default/viewer.htm#ormpug_decomp_examples02.htm, an example made by SAS
# corporation. 24 tasks assigned to 8 agents.
# =============================================================================
cost = t(matrix(c(
  8,18,22,5,11,11,22,11,17,22,11,20,13,13,7,22,15,22,24,8,8,24,18,8,24,14,11,
  15,24,8,10,15,19,25,6,13,10,25,19,24,13,12,5,18,10,24,8,5,22,22,21,22,13,
  16,21,5,25,13,12,9,24,6,22,24,11,21,11,14,12,10,20,6,13,8,19,12,19,18,10,21,
  5,9,11,9,22,8,12,13,9,25,19,24,22,6,19,14,25,16,13,5,11,8,7,8,25,20,24,20,11,
  6,10,10,6,22,10,10,13,21,5,19,19,19,5,11,22,24,18,11,6,13,24,24,22,6,22,5,14,
  6,16,11,6,8,18,10,24,10,9,10,6,15,7,13,20,8,7,9,24,9,21,9,11,19,10,5,23,20,5,
  21,6,9,9,5,12,10,16,15,19,18,20,18,16,21,11,12,22,16,21,25,7,14,16,10),
  ncol = 8))
profit = t(matrix(c(
  25,23,20,16,19,22,20,16,15,22,15,21,20,23,20,22,19,25,25,24,21,17,23,17,16,
  19,22,22,19,23,17,24,15,24,18,19,20,24,25,25,19,24,18,21,16,25,15,20,20,18,
  23,23,23,17,19,16,24,24,17,23,19,22,23,25,23,18,19,24,20,17,23,23,16,16,15,23,
  15,15,25,22,17,20,19,16,17,17,20,17,17,18,16,18,15,25,22,17,17,23,21,20,24,22,
  25,17,22,20,16,22,21,23,24,15,22,25,18,19,19,17,22,23,24,21,23,17,21,19,19,17,
  18,24,15,15,17,18,15,24,19,21,23,24,17,20,16,21,18,21,22,23,22,15,18,15,21,22,
  15,23,21,25,25,23,20,16,25,17,15,15,18,16,19,24,18,17,21,18,24,25,18,23,21,15,
  24,23,18,18,23,23,16,20,20,19,25,21), ncol = 8))
budget = c(36, 35, 38, 34, 32, 34, 31, 34)


# Intel CPU i7-4770 3.4GHz, g++ '-Ofast', 64-bit Windows 7:
system.time({sol = FLSSS::auxGAPbb(
  cost, profit, budget, maxCore = 2, tlimit = 4, ub = "MT",
  greedyBranching = FALSE, optim = "max")})
# user  system elapsed
# 0.02    0.00    0.01
# The elapsed time is about 1\% of that reported by the SAS proc with 8
# threads, although its hardware configuration is unknown.


system.time({sol2 = FLSSS::auxGAPbb(
  cost, profit, budget, maxCore = 2, tlimit = 4, ub = "MT",
  greedyBranching = TRUE, optim = "max")})
sol[c("nodes", "bkpSolved")]  # 4526, 14671, can be different.
sol2[c("nodes", "bkpSolved")] # 4517, 13115, can be different.
# Greedy branching may lower the numbers of branching nodes and
# knapsack problems to solve.



\donttest{
# =============================================================================
# Play random numbers.
# =============================================================================
set.seed(22) # A nontrivial instance searched via changing random seeds.
             # RNG in R 3.5.1 for Windows.
Nagent = 20L; Ntask = 200L
cost = matrix(runif(Nagent * Ntask, 1e3, 1e6), nrow = Nagent)
profit = matrix(abs(rnorm(Nagent * Ntask, 1e6, 1e6)) + 1000, nrow = Nagent)
budget = apply(cost, 1, function(x) runif(1, min(x), sum(x) / 2))


# Intel CPU i7-4770 3.4GHz, g++ '-Ofast', 64-bit Windows 7.
system.time({sol1 = FLSSS::auxGAPbb(
  cost, profit, budget,
  maxCore = 1, multhreadOn = "KPs",
  tlimit = 3600, ub = "MT", greedyBranching = TRUE, optim = "max")})
# user  system elapsed
# 9.17    0.00    9.19


# Multithread knapsack problems at each branching node.
# This does not accelerate the speed at all because threading overheads
# are overwhelming.
system.time({sol2 = FLSSS::auxGAPbb(
  cost, profit, budget,
  maxCore = 7, multhreadOn = "KPs",
  tlimit = 3600, ub = "MT", greedyBranching = TRUE, optim = "max")})
# user   system elapsed
# 39.02    5.24   11.12


# Multithread nodes.
system.time({sol3 = FLSSS::auxGAPbb(
  cost, profit, budget,
  maxCore = 7, multhreadOn = "nodes", threadLoad = 32L,
  tlimit = 3600, ub = "MT", greedyBranching = TRUE, optim = "max")})
#  user  system elapsed
# 14.62    0.00    2.13
}
}

 