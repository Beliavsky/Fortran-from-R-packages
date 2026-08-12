\name{subsetsum}
\alias{subsetsum}
\alias{sss_test}
\title{
  Subset Sum Problem
}
\description{
  Subset sum routine for positive integers.
}
\usage{
subsetsum(S, t, method = "greedy")

sss_test(S, t)
}
\arguments{
  \item{S}{vector of positive integers.}
  \item{t}{target value, bigger than all items in \code{S}.}
  \item{method}{can be ``greedy'' or ``dynamic'', where ``dynamic'' stands
                for the dynamic programming approach.}
}
\details{
  \code{subsetsum} is searching for a set of elements in \code{S} that
  sum up to \code{t} by continuously adding more elements of \code{S}.

  It is not required that \code{S} is decreasingly sorted. But for reasons
  of efficiency and smaller execution times it is urgently recommended to
  sort the item set in decreasing order. See the examples to find out how
  to handle your data.

  The first components will be preferred, i.e., if \code{S} is decreasing,
  the sum with larger elements will be found, if increasing, the sum with
  smaller elements. Because of timing considerations, the default is to
  sort decreasingly before processing.

  The dynamic method may be faster for large sets, but will also require
  much more memory if the target value is large. 

  \code{sss_test} will find the biggest number below or equal to \code{t}
  that can be expressed as a sum of items in \code{S}. It will not return
  any indices. It can be quite fast, though it preprocesses the set \code{S}
  to be sorted decreasingly, too. 
}
\value{
  List with the target value, if reached, and vector of indices of elements
  in \code{S} that sum up to \code{t}.

  If no solution is found, the dynamic method will return indices for the
  largest value below the target, the greedy method witll return NULL.

  \code{sss_test} will simply return maximum sum value found.
}
\note{
  A compiled version -- and much faster, in Fortran -- can be found in
  package 'knapsack' (R-Forge, project 'optimist') as \code{subsetsum}.
  A recursive version, returning *all* solutions, is much too slow in R,
  but is possible in Julia and can be asked from the author.
}
\author{
  HwB  email: <hwborchers@googlemail.com>
}
\references{
  Horowitz, E., and S. Sahni (1978). Fundamentals of Computer Algorithms.
  Computer Science Press, Rockville, ML.
}
\seealso{
  \code{\link{maxsub}}
}
\examples{
t <- 5842
S <- c(267, 493, 869, 961, 1000, 1153, 1246, 1598, 1766, 1922)

# S is not decreasingly sorted, so ...
o  <- order(S, decreasing = TRUE)
So <- S[o]                          # So is decreasingly sorted

sol <- subsetsum(So, t)             # $inds:  2 4 6 7 8  w.r.t.  So
is  <- o[sol$inds]                  # is:     9 7 5 4 3  w.r.t.  S
sum(S[is])                          # 5842

\dontrun{
amount <- 4748652
products <- 
c(30500,30500,30500,30500,42000,42000,42000,42000,
  42000,42000,42000,42000,42000,42000,71040,90900,
  76950,35100,71190,53730,456000,70740,70740,533600,
  83800,59500,27465,28000,28000,28000,28000,28000,
  26140,49600,77000,123289,27000,27000,27000,27000,
  27000,27000,80000,33000,33000,55000,77382,48048,
  51186,40000,35000,21716,63051,15025,15025,15025,
  15025,800000,1110000,59700,25908,829350,1198000,1031655)

# prepare set
prods <- products[products <= amount]  # no elements > amount
prods <- sort(prods, decreasing=TRUE)  # decreasing order

# now find one solution
system.time(is <- subsetsum(prods, amount))
#  user  system elapsed 
# 0.030   0.000   0.029 

prods[is]
#  [1]   70740   70740   71190   76950   77382   80000   83800
#  [8]   90900  456000  533600  829350 1110000 1198000

sum(prods[is]) == amount
# [1] TRUE

# Timings:
#             unsorted   decr.sorted
# "greedy"      22.930         0.030    (therefore the default settings)
# "dynamic"      2.515         0.860    (overhead for smaller sets)
# sss_test       8.450         0.040    (no indices returned)
}
}
\keyword{ discrete-optimization }
