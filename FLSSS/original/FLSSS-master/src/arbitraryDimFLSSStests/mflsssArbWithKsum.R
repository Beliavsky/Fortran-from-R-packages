


Rcpp::sourceCpp("src/mflsssArb.cpp", verbose = T)
Rcpp::sourceCpp("src/arbitraryDimFLSSStests/findBound.cpp", verbose = T)
source("src/arbitraryDimFLSSStests/stringMatTo64bitIntMat.R")
source("R/numStrings.R")




set.seed(123456)
options(scipen = 999)


# Rcpp::sourceCpp("src/mflsssArb.cpp", verbose = T)
iter = 0L
while (T) {


  iter = iter + 1L
  cat(iter, "")


  d = sample(5, 1) + 1L
  N = sample(1:50, 1)
  len = min(N, sample(1:N, 1))
  if (runif(1) < 0.5) len = max(N - len, 1L)
  roundN = sample(5, 1)


  X = matrix(round(runif(d * N, -1, 1), roundN), nrow = N)
  if (runif(1) < 0.5) sol = sort(sample(N, len)) else sol = integer(0)
  save(N, d, len, X, sol, roundN, file = "debug.Rdata")


  # load("debug.Rdata")
  if (length(sol) == 0)
  {
    tmp = colSums(X[sample(N, len), , drop = F])
    tmp2 = colSums(X[sample(N, len), , drop = F])
    target = unlist(mapply(function(x, y)
    {
      z = c(x, y)
      round(runif(1, min(z), max(z)), roundN)
    }, tmp, tmp2))
  } else {
    target = round(colSums(X[sol, , drop = F]), roundN)
  }


  Xc = matrix(as.character(X), nrow = N)
  targetc = as.character(target)


  tmpX = round(X * 10 ^ roundN)
  tmpXmin = apply(tmpX, 2, function(x) min(x))
  tmpX = t(t(tmpX) - tmpXmin)
  tmpTarget = round(target * 10 ^ roundN)
  tmpTarget = round(tmpTarget - tmpXmin * len)
  tmpXrank = apply(tmpX, 2, function(x) rank(x, ties.method = "first"))
  colOdr = order(colSums(t(tmpXrank) %*% tmpXrank))
  Xc = Xc[, colOdr, drop = F]
  targetc = targetc[colOdr]
  tmpX = tmpX[, colOdr, drop = F]
  tmpTarget = as.integer(tmpTarget[colOdr])
  rowOdr = do.call(order, as.data.frame(tmpX[, ncol(tmpX):1, drop = F]))
  tmpX = matrix(as.integer(tmpX[rowOdr, , drop = F]), nrow = nrow(tmpX))
  Xc = Xc[rowOdr, , drop = F]
  tmpsol = sort(match(sol, rowOdr))


  ksumK = sample(3:6, 1)
  choice1 = runif(1) < 0.5
  if (choice1) tmphashtarget = NULL else tmphashtarget = targetc


  choice2 = runif(1) < 0.5
  if (choice2) ksumtablesupply = NULL else {
    ksumtablesupply = ksumHash(
    ksumK = ksumK, V = Xc, ksumTableSizeScaler = 30,
    target = tmphashtarget, len = len, approxNinstance = 1000,
    verbose = F, maxCore = 7)
  }


  system.time({arbflssssol = arbFLSSS(
    len = len, V = Xc, target = targetc, givenKsumTable = ksumtablesupply,
    solutionNeed = 1e9, maxCore = 10, approxNinstance = 1000,
    ksumK = ksumK, tlimit = 1e9, ksumTableSizeScaler = 30L, verbose = F)})


  choice3 = runif(1) < 0.5
  if (choice3) decomposeArbFLSSSksumTable = ksumtablesupply else
    decomposeArbFLSSSksumTable = NULL


  tmp = decomposeArbFLSSS(
    len = len, V = Xc, target = targetc, approxNinstance = 1000,
    ksumTable = decomposeArbFLSSSksumTable, ksum = ksumK,
    ksumTableSizeScaler = 30L, verbose = F, maxCore = 10)
  existsol = tmp$solutionsFound
  minedsol = unlist(lapply(tmp$arbFLSSSobjects, function(x)
  {
    arbFLSSSobjRun(
      x, solutionNeed = 1e9, tlimit = 1e9, ksum = ksumK,
      ksumTableSizeScaler = 30, verbose = F)
  }), recursive = F)


  decompsol = sort(unlist(lapply(c(existsol, minedsol), function(x)
  {
    paste0(sort(x), collapse = "-")
  })))
  arbflssssol = sort(unlist(lapply(arbflssssol, function(x)
  {
    paste0(sort(x), collapse = "-")
  })))


  arbflsssAndDecompEqual = all(decompsol == arbflssssol)
  if (!arbflsssAndDecompEqual)
  {
    cat("arbflsss and decomp not equal!\n"); break
  }


  correctRst = FLSSS::mFLSSSpar(
    maxCore = 10L, len = len, mV = tmpX, mTarget = tmpTarget,
    mME = numeric(ncol(tmpX)) + 0.3, solutionNeed = 1e9L,
    tlimit = 1e9, useBiSrchInFB = F, avgThreadLoad = 100L)


  correctRstSol = sort(unlist(lapply(correctRst, function(x)
  {
    paste0(sort(x), collapse = "-")
  })))


  arbflsssAndCorrectEqual = all(correctRstSol == arbflssssol)
  if (!arbflsssAndCorrectEqual)
  {
    cat("arbflsss And Correct not Equal!\n"); break
  }


}



































