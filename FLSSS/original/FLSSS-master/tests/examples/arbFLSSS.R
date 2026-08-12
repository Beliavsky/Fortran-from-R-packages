# ==============================================================================
# Mine subsets with arbitrary precision and magnitude.
# ==============================================================================
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


system.time({
  rst = FLSSS::arbFLSSS(
    len, V = as.matrix(V), target, solutionNeed = 1, maxCore = 2,
    tlimit = 10, ksumK = 0, verbose = TRUE)
})


# Validation.
all(unlist(lapply(rst, function(x) FLSSS:::addNumStrings(V[x]))) == target)




# ==============================================================================
# Mine in a multidimensional set.
# ==============================================================================
set.seed(2)
d = 4L # Set dimension.
N = 50L # Set size.
len = 10L # Subset size.
roundN = 4L # For rounding the numeric values before conversion to strings.


V = matrix(round(runif(N * d, -1, 1), roundN), nrow = N) # Make superset.
options(scipen = 999) # Ensure numeric-to-string conversion does not
# produce strings like "2e-3".
Vstr = matrix(as.character(V), nrow = N)


sol = sample(N, len) # Make a solution.
target = round(colSums(V[sol, ]), roundN) # Target subset sum.
targetStr = as.character(target)


system.time({
  rst = FLSSS::arbFLSSS(
    len = len, V = Vstr, target, givenKsumTable = NULL, tlimit = 60,
    solutionNeed = 1e9, maxCore = 2, ksumK = 4, verbose = TRUE)
})


# Validation.
all(unlist(lapply(rst, function(x)
{
  apply(Vstr, 2, function(u) FLSSS:::addNumStrings(u[x]))
})) == targetStr)




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




































