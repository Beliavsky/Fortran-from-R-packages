


set.seed(42)
d = 5L # Set dimension.
N = 30L # Set size.
len = 10L # Subset size.
roundN = 2L # For rounding the numeric values before conversion to strings.


V = matrix(round(runif(N * d, -1e5, 1e5), roundN), nrow = N) # Make superset.
sol = sample(N, len) # Make a solution.
target = round(colSums(V[sol, ]), roundN) # Target subset sum.


options(scipen = 999) # Ensure numeric => string conversion does not
# produce strings like 2e-3.
Vstr = matrix(as.character(V), nrow = N) # String version of V.
targetStr = as.character(target)


system.time({
  theDecomposed = FLSSS::decomposeArbFLSSS(
    len = len, V = Vstr, target = targetStr, approxNinstance = 1000,
    maxCore = 2, ksumTable = NULL, ksumK = 4, verbose = TRUE)
})


# Run the objects sequentially.
rst = unlist(lapply(theDecomposed$arbFLSSSobjects, function(x)
{
  FLSSS::arbFLSSSobjRun(x, solutionNeed = 1e9, tlimit = 5, verbose = FALSE)
}), recursive = FALSE)
str(rst)


















