\name{FLSSS}
\alias{FLSSS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
One-dimensional Subset Sum given error threshold
}
\description{
Given subset size \code{len}, sorted superset \code{v}, subset sum \code{target} and error \code{ME}, find at least \code{solutionNeed} index (integer) vector(s) \code{x}, such that \code{target - ME <= sum(v[x]) <= target + ME}. To mine subsets that sum in a given range, set \code{target} to the midpoint and \code{ME} to half of the range width.
}
\usage{
FLSSS(
  len,
  v,
  target,
  ME,
  solutionNeed = 1L,
  LB = 1L : len,
  UB = (length(v) - len + 1L) : length(v),
  viaConjugate = FALSE,
  tlimit = 60,
  useBiSrchInFB = FALSE,
  NfractionDigits = Inf
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{len}{
An integer as the subset size: \code{0 <= len < length(v)}. If \code{len == 0}, \code{FLSSS()} mines subets without size restriction. \code{len <- 0} would be most likely slower than looping \code{len} over \code{1 : (length(v) - 1)}. See Details.
}
  \item{v}{
A sorted numeric vector, the superset. \code{v} can be negative and nonunique.
}
  \item{target}{
A numeric value, the subset sum target.
}
  \item{ME}{
A positive numeric value, the error threshold.
}
  \item{solutionNeed}{
An integer, the least number of solutions wanted. If the function returns fewer solutions, either \code{tlimit} is up or less than \code{solutionNeed} solutions exist. The function may also return more than \code{solutionNeed} solutions.
}
  \item{LB}{
An integer vector of size \code{len} as the lower bounds of the solution space: for any solution \code{x}, \code{LB[i] <= x[i]}. Custom \code{LB} should be no less than \code{1L : len} element-wisely. Every element in \code{v} should be within the range enclosed by \code{LB} and \code{UB}.
}
  \item{UB}{
An integer vector of size \code{len} as the upper bounds of the solution space: for any solution \code{x}, \code{x[i] <= UB[i]}. Custom \code{UB} should be no greater than \code{(length(v)} \code{-} \code{len} \code{+} \code{1L} \code{)} \code{:} \code{length(v)} element-wisely. Every element in \code{v} should be within the range enclosed by \code{LB} and \code{UB}.
}
  \item{viaConjugate}{
A boolean value. If \code{TRUE}, \code{FLSSS()} mines susbets of size \code{length(v) - len} that sum to \code{sum(v)} \code{- target} with the same \code{ME}. Let \code{x} be the integer vector indexing a qualified subset. \code{FLSSS()} returns \code{(1L : length(v))[-x]}. Simulations show that \code{FLSSS()} often finds the first qualified conjugate subset faster if \code{len} is much less than \code{length(v) / 2}.
}
  \item{tlimit}{
A numeric value. Enforce function to return in \code{tlimit} seconds.
}
  \item{useBiSrchInFB}{
A boolean value. If \code{TRUE}, the function performs binary search for index bounds in the auxiliary triangle matrix of continuous sequence sums. This argument is mainly for research. Simulations show binary search has no major advantage over linear search due to caching mechanisms. The advantage may be pronounced if \code{length(v)} is substantial ( > 10000) while \code{len} is small ( < 5).
}
  \item{NfractionDigits}{
An integer, the maximum number of fractional digits of all elements in \code{v}. Internally, \code{v}, \code{target} and \code{ME} are multiplied by \code{10 ^ NfractionDigits}, and then converted as integer values before mining. The default \code{Inf} prevents such conversion. The goal is eliminate
}
}
\details{
If \code{len == 0}, \code{FLSSS()} would (1) reset \code{len} to \code{length(v)}, (2) pad \code{len} zeros at the beginning of \code{v} and sort \code{v}, (3) search for size-\code{len} subsets, and (4) for an index vector that represents a subset, erases elements pointing to zeros in \code{v}. See the \href{https://arxiv.org/abs/1612.04484}{package documentation} for more details.
}
\value{
A list of index vectors.
}
\examples{
# =====================================================================================
# Example I: play random numbers.
# =====================================================================================
# rm(list = ls()); gc()
subsetSize = 200L
supersetSize = 1000L
superset = 10000 * sort(rnorm(supersetSize) ^ 3 + 2 * runif(supersetSize) ^ 2 +
           3 * rgamma(supersetSize, 5, 1) + 4)
subsetSum = runif(1, sum(superset[1L : subsetSize]), sum(superset[(supersetSize -
            subsetSize + 1L) : supersetSize]))
subsetSumError = 1e-3


# Mine 3 subsets
rst1 = FLSSS::FLSSS(len = subsetSize, v = superset, target = subsetSum,
                    ME = subsetSumError, solutionNeed = 3, tlimit = 4)


# Mine 3 subsets via solving the conjugate problem
rst2 = FLSSS::FLSSS(len = subsetSize, v = superset, target = subsetSum,
                    ME = subsetSumError, solutionNeed = 3, tlimit = 4,
                    viaConjugate = TRUE)


# Verify uniqueness
cat("rst1 number of solutions =",
    length(unique(lapply(rst1, function(x) sort(x)))), "\n")
cat("rst2 number of solutions =",
    length(unique(lapply(rst2, function(x) sort(x)))), "\n")


# Verify solutions
if(length(rst1) > 0)
  all(unlist(lapply(rst1, function(x)
    abs(sum(superset[x]) - subsetSum) <= subsetSumError)))
if(length(rst2) > 0)
  all(unlist(lapply(rst2, function(x)
    abs(sum(superset[x]) - subsetSum) <= subsetSumError)))


# Mine 3 subsets in bounded solution space.
# Make up the lower and upper bounds for the solution space:
tmp = sort(sample(1L : supersetSize, subsetSize))
tmp2 = sort(sample(1L : supersetSize, subsetSize))
lowerBounds = pmin(tmp, tmp2)
upperBounds = pmax(tmp, tmp2)
rm(tmp, tmp2)


# 'FLSSS()' does not work if there are elements not under the hood of
# lowerBounds + upperBounds. Exclude those elements:
remainIndex = unique(unlist(apply(cbind(lowerBounds, upperBounds), 1,
  function(x) x[1] : x[2])))
lowerBounds = match(lowerBounds, remainIndex)
upperBounds = match(upperBounds, remainIndex)
superset = superset[remainIndex]


# Plant a subset sum:
solution = integer(subsetSize)
solution[1] = sample(lowerBounds[1] : upperBounds[1], 1)
for(i in 2L : subsetSize)
{
  l = max(lowerBounds[i], solution[i - 1] + 1L)
  u = upperBounds[i]
  if(l == u) solution[i] = u
  else solution[i] = sample(l : u, 1)
}
subsetSum = sum(superset[solution])
subsetSumError = abs(subsetSum) * 0.01 # relative error within 1%
rm(solution)


rst3 = FLSSS::FLSSS(len = subsetSize, v = superset, target = subsetSum,
                    ME = subsetSumError, solutionNeed = 2, tlimit = 4,
                    LB = lowerBounds, UB = upperBounds, viaConjugate = TRUE)


print(length(rst3))


# Verify solutions
if(length(rst3) > 0)
  cat(all(unlist(lapply(rst3, function(x)
    abs(sum(superset[x]) - subsetSum) <= subsetSumError))), "\n")




# =====================================================================================
# Example II: mine a real-world dataset.
# =====================================================================================
# rm(list = ls()); gc()
superset = c(
  -1119924501, -793412295, -496234747,  -213654767,   16818148,   26267601,   26557292,
     27340260,   28343800,   32036573,    32847411,   34570996,   34574989,   43633028,
     44003100,   47724096,   51905122,    52691025,   53600924,   56874435,   58207678,
     60225777,   60639161,   60888288,    60890325,   61742932,   63780621,   63786876,
     65167464,   66224357,   67198760,    69366452,   71163068,   72338751,   72960793,
     73197629,   76148392,   77779087,    78308432,   81196763,   82741805,   85315243,
     86446883,   87820032,   89819002,    90604146,   93761290,   97920291,   98315039,
    310120088, -441403864, -548143111,  -645883459, -149110919,  305170449, -248934805,
  -1108320430, -527806318, -192539936, -1005074405, -101557770, -156782742, -285384687,
   -418917176,   80346546, -273215446,  -552291568,   86824498,  -95392618, -707778486)
superset = sort(superset)
subsetSum = 139254953
subsetSumError = 0.1


# Find a subset of size 10.
subsetSize = 10L
rst = FLSSS::FLSSS(len = subsetSize, v = superset, target = subsetSum,
                   ME = subsetSumError, solutionNeed = 1, tlimit = 4)
# Verify:
all(unlist(lapply(rst, function(x)
  abs(sum(superset[x]) - subsetSum) <= subsetSumError)))


# Find a subset without size specification.
rst = FLSSS::FLSSS(len = 0, v = superset, target = subsetSum,
                   ME = subsetSumError, solutionNeed = 1, tlimit = 4)
# Verify:
all(unlist(lapply(rst, function(x)
  abs(sum(superset[x]) - subsetSum) <= subsetSumError)))


# Find a subset via looping subset size over 2L : (length(v)).
for(len in 2L : length(superset))
{
  rst = FLSSS::FLSSS(len = subsetSize, v = superset, target = subsetSum,
                     ME = subsetSumError, solutionNeed = 1, tlimit = 4)
  if(length(rst) > 0) break
}
# Verify:
all(unlist(lapply(rst, function(x)
  abs(sum(superset[x]) - subsetSum) <= subsetSumError)))


# Find as many qualified susbets as possible in 2 seconds
rst = FLSSS::FLSSS(len = subsetSize, v = superset, target = subsetSum,
                   ME = subsetSumError, solutionNeed = 999999L, tlimit = 2)
cat("Number of solutions =", length(rst), "\n")


# Verify:
all(unlist(lapply(rst, function(x)
  abs(sum(superset[x]) - subsetSum) <= subsetSumError)))




# =====================================================================================
# Example III: solve a special knapsack problem.
# Given the knapsack's capacity, the number of catagories, the number of items in each
# catagory, select the least number of items to fulfill at least 95% of the knapsack's
# capacity.
# =====================================================================================
# rm(list = ls()); gc()
capacity = 361
catagories = LETTERS[1L : 10L] # A, B, ..., J, 10 catagories
catagoryMasses = round(runif(length(catagories)) * 20 + 1)
catagoryItems = sample(1L : 20L, length(catagories))


itemLabel = unlist(mapply(function(x, i) rep(i, x), catagoryItems, catagories))
itemMasses = unlist(mapply(function(x, i) rep(x, i), catagoryMasses, catagoryItems))
vorder = order(itemMasses)
itemLabel = itemLabel[vorder]


superset = itemMasses[vorder]
rate = 0.95
subsetSum = (capacity * rate + capacity) / 2
subsetSumError = capacity - subsetSum
for(subsetSize in 1L : length(itemMasses))
{
  rst = FLSSS::FLSSS(len = subsetSize, v = superset, target = subsetSum,
                     ME = subsetSumError, solutionNeed = 1, tlimit = 4)
  if(length(rst) > 0) break
}


# There may exist no qualified subsets. One can lower 'rate' until a solution
# shows up.
if(length(rst) == 0L)
{
  cat("No solutions. Please lower rate and rerun.\n")
} else
{
  cat("A solution:\n")
  print(table(itemLabel[rst[[1]]]))
}


# rm(list = ls()); gc()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
 