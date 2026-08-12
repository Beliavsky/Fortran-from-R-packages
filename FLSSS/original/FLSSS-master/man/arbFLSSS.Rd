\name{arbFLSSS}
\alias{arbFLSSS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Multidimensional exact subset sum in arbitrary precision and magnitude
}
\description{
Given a multidimensional set and a subset size, find one or more subsets whose elements sum up to a given target.
}
\usage{
arbFLSSS(
  len,
  V,
  target,
  givenKsumTable,
  solutionNeed = 1L,
  maxCore = 7L,
  tlimit = 60,
  approxNinstance = 1000L,
  ksumK = 4L,
  ksumTableSizeScaler = 30L,
  verbose = TRUE
  )
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{len}{
An integer as the subset size. \code{1 <= len <= nrow(V)}.
}
  \item{V}{
A string matrix as the superset. Rows are elements.
}
  \item{target}{
A string vector as the target subset sum. \code{length(target) == ncol(V)}.
}
  \item{givenKsumTable}{
Either \code{NULL} or the return value from \code{ksumHash()}. See argument \code{ksumK} for the preliminaries. If \code{NULL}, the function will compute and hash k-sums depending on \code{ksumK} before mining the subsets. Otherwise it will use \code{givenKsumTable} as the lookup table and ignore arguments \code{ksumK} and \code{ksumTableSizeScaler}.
}
  \item{solutionNeed}{
An integer. How many solutions are wanted. Default \code{solutionNeed = 1}.
}
  \item{maxCore}{
An integer as the maximum threads to invoke. Better not exceed the number of logical processors on the platform. Default \code{maxCore = 7}.
}
  \item{tlimit}{
A numeric value as the time limit (seconds). Default \code{tlimit = 60}.
}
  \item{approxNinstance}{
An integer. The problem will be decomposed into about \code{approxNinstance} subproblems solved independently by the threads. Default \code{approxNinstance = 1000}. \code{approxNinstance} is better to be much higher than argument \code{maxCore} since time costs of the subproblems are unknown and probably vary greatly.
}
  \item{ksumK}{
An integer. If \code{ksumK < 3}, no k-sum accelerator will be built. For example, if \code{ksumK = 5}, then the sums of all combinations of 3 elements (3-sums), the sums of all combinations of 4 elements (4-sums), and the sums of all combinations of 5 elements (5-sums) in \code{V}, are pre-computed and hashed into a Bloom filter variant. This filter thus contains 3 lookup tables responding to the 3-sums, 4-sums and 5-sums. During the main course of mining, if any set is reduced to one of size 3, 4, or 5, the set's associated target sum will be hashed and looked up in the filter. Not existing would imply the target sum is unreachable and the set can be discredited immediately. This typically generates massive speedup. For \code{ksumK < 3}, such filtering is not meaningful and thus not performed. A high \code{ksumK} coupled with a large superset \code{V} however is prone to memory overflow or extremely time consuming. \code{ksumK} will be upper-bounded by subset size \code{len} internally. Default \code{ksumK = 4}.
}
  \item{ksumTableSizeScaler}{
An integer for determining size of the k-sum lookup table in the filter described above. For example, a set of size 21 has 1330 3-element subsets. If \code{ksumTableSizeScaler = 10}, then around 13300 bits will be allocated for the 3-sum lookup table. The exact number of bits is 14033 + 7, where 14033 is the lowest element greater than 13300 in a prime array defined in GCC's STL of hashing policy, and 7 is to make up the last byte. Default  \code{ksumTableSizeScaler = 30}. Higher \code{ksumTableSizeScaler} means lower chance of hash collision, thus higher efficiency.
}
  \item{verbose}{
A boolean value. \code{TRUE} prints the computing progress. Default \code{TRUE}.
}
}
\details{
New users might want to check out \code{FLSSS()} or \code{mFLSSSpar()} first.

String matrix \code{V} is maximally compressed into an integer set of size \code{nrow(V)}. Dimensionality of the set will be printed given \code{verbose = TRUE}. Each set element is a huge integer comprising many 64-bit buffers. Addition and subtraction of the huge integers call \code{mpn_add_n()} and \code{mpn_sub_n()} from the GNU Multiple Precision Arithmetic Library (GMP) if the system has it, otherwise they are performed by customized algorithms.

After the initial problem is decomposed, the smaller problems can collectively offer a pair of index lower and upper bounds. The k-subsets outside the bounds are not necessarily considered for building the k-sum accelerator.

See comparisons between this function and \code{FLSSS()}, \code{mFLSSSpar()} in Examples.
}
\value{
A list of index vectors as solutions.
}


\examples{
set.seed(1)
N = 200L # Superset size.
len = 20L # Subset size.
V = sapply(1:N, function(i) # Generate a set where every "number" has at most
  # 100 digits.
{
  a = 0:9
  left = sample(a, size = sample(50, 1), replace = TRUE)
  right = sample(a, size = sample(50, 1), replace = TRUE)
  x = paste0(paste0(left, collapse = ""), ".", paste0(right, collapse = ""))
  if (runif(1) < 0.5) x = paste0("-", x) # Randomly add a negative sign.
  x
})
str(V)


sol = sample(N, len) # Make a solution.
target = FLSSS::addNumStrings(V[sol]) # An unexposed helper function.


# Uncomment the following block for testing the code since CRAN check uses 
# Valgrind, which takes over 5s to run.
# system.time({
#   rst = FLSSS::arbFLSSS(
#     len = len, V = Vstr, target, givenKsumTable = NULL, tlimit = 60,
#     solutionNeed = 1e9, maxCore = 2, ksumK = 4, verbose = TRUE)
# })
# 
# 
# # Validation.
# all(unlist(lapply(rst, function(x)
# {
#   apply(Vstr, 2, function(u) FLSSS:::addNumStrings(u[x]))
# })) == targetStr)



# ==============================================================================
# Mine in a multidimensional set.
# ==============================================================================
set.seed(2)
d = 4L # Set dimension.
N = 50L # Set size.
len = 10L # Subset size.
roundN = 4L # For rounding the numeric values before conversion to strings.


V = matrix(round(runif(N * d, -1, 1), roundN), nrow = N) # Make superset.
optionSave = options()
options(scipen = 999) # Ensure numeric-to-string conversion does not
# produce strings like "2e-3".
Vstr = matrix(as.character(V), nrow = N)


sol = sample(N, len) # Make a solution.
target = round(colSums(V[sol, ]), roundN) # Target subset sum.
targetStr = as.character(target)


# Uncomment the following block for testing the code since CRAN check uses 
# Valgrind, which takes over 5s to run.
# system.time({
#   rst = FLSSS::arbFLSSS(
#     len = len, V = Vstr, target, givenKsumTable = NULL, tlimit = 60,
#     solutionNeed = 1e9, maxCore = 2, ksumK = 4, verbose = TRUE)
# })
# 
# 
# # Validation.
# all(unlist(lapply(rst, function(x)
# {
#   apply(Vstr, 2, function(u) FLSSS:::addNumStrings(u[x]))
# })) == targetStr)




# # ============================================================================
# # Compare arbFLSSS() and FLSSS(). Example takes more than 2 seconds. The 
# # section has some analysis of the algorithms.
# # ============================================================================
# set.seed(3)
# N = 100L # Superset size.
# len = 20L # Subset size.
# roundN = 5L # For rounding the numeric values.
# V = sort(round(100000 * runif(N, -1, 1), roundN)) # Create superset.
# sol = sort(sample(N, len)) # Make a solution.
# target = round(sum(V[sol]), roundN)
# error = 3e-6 # Effectively demands the target sum to be exactly matched
# # since roundN = 5.
# 
# 
# system.time({
#   FLSSSrst = FLSSS::FLSSS(
#     len, V, target, ME = error, solutionNeed = 2, tlimit = 60)
# })
# # It may seem counter-intuitive that this takes much longer than the instance
# # with N = 1000 and len = 200L --- the 1st example in the help page of
# # FLSSS(). Note the time cost is closely related to the "rarity" of
# # solutions. A larger superset or subset could mean more element combinations
# # that can sum into the given range, thus more solutions and easier to mine.
# 
# 
# # Validate the results.
# all(abs(unlist(lapply(FLSSSrst, function(x) sum(V[x]))) - target) <= error)
# 
# 
# options(scipen = 999)
# Vstr = as.matrix(as.character(V))
# targetStr = as.character(target)
# # Use 1 thread for a fair comparison with FLSSS() since the latter is
# # single-threaded. Use no k-sum accelerator.
# system.time({
#   arbFLSSSrst = FLSSS::arbFLSSS(
#     len, V = Vstr, target = targetStr, solutionNeed = 2, maxCore = 1,
#     ksumK = 0, verbose = TRUE, approxNinstance = 1000, tlimit = 60)
# })
# # Timing is higher than FLSSS() because arbFLSSS()'s objective
# # is not just solving unidimensional problem.
# 
# 
# # Validation.
# all(abs(unlist(lapply(arbFLSSSrst, function(x) sum(V[x]))) - target) <= error)
# 
# 
# # Use 4-sum accelerator. Massive speedup.
# system.time({
#   arbFLSSSrst = FLSSS::arbFLSSS(
#     len, Vstr, targetStr, solutionNeed = 2, maxCore = 1, ksumK = 4,
#     verbose = FALSE, approxNinstance = 1000, tlimit = 60)
# })
# 
# 
# # Validation.
# all(abs(unlist(lapply(arbFLSSSrst, function(x) sum(V[x]))) - target) <= error)




# # ============================================================================
# # Compare arbFLSSS() and mFLSSSpar(). Example takes more than 2 seconds. The
# # section contains some analysis of the algorithms.
# # ============================================================================
# set.seed(4)
# d = 5L # Set dimension.
# N = 60L # Set size.
# len = 10L # Subset size.
# roundN = 2L # For rounding the numeric values before conversion to strings.
# 
# 
# V = matrix(round(runif(N * d, -1e5, 1e5), roundN), nrow = N) # Make superset.
# sol = sample(N, len) # Make a solution.
# target = round(colSums(V[sol, ]), roundN) # Target subset sum.
# error = rep(2e-3, d) # Effectively demands the target sum to be exactly
# # matched since roundN = 2.
# 
# 
# system.time({
#   mFLSSSparRst = FLSSS::mFLSSSpar(
#     maxCore = 7, len = len, mV = V, mTarget = target, mME = error,
#     avgThreadLoad = 20, solutionNeed = 1, tlimit = 60)
# })
# 
# 
# # Validation.
# all(unlist(lapply(mFLSSSparRst, function(x)
#   abs(colSums(V[x, , drop = FALSE]) - target) <= error)))
# 
# 
# options(scipen = 999) # Ensure numeric => string conversion does not
# # produce strings like 2e-3.
# Vstr = matrix(as.character(V), nrow = N) # String version of V.
# targetStr = as.character(target)
# 
# 
# # Use no k-sum accelerator.
# system.time({
#   arbFLSSSrst = FLSSS::arbFLSSS(
#     len = len, V = Vstr, target = targetStr, givenKsumTable = NULL,
#     tlimit = 60, solutionNeed = 1, maxCore = 7, ksumK = 0, verbose = TRUE)
# })
# 
# 
# # Validation.
# all(unlist(lapply(arbFLSSSrst, function(x)
#   abs(colSums(V[x, , drop = FALSE]) - target) <= error)))
# 
# 
# # Use 5-sum accelerator. Massive speedup.
# system.time({
#   arbFLSSSrst = FLSSS::arbFLSSS(
#     len = len, V = Vstr, target = targetStr, givenKsumTable = NULL,
#     tlimit = 60, solutionNeed = 1, maxCore = 7, ksumK = 5, verbose = TRUE)
# })


options(optionSave)
}















