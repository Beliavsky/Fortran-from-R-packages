

source("src/arbitraryDimFLSSStests/findBound.R")


mflsssArbR = function(len, v, target, lb, ub, hopeV)
{
  rst = list()


  if (len == 1L)
  {
    b = lb[1]:ub[1]
    ind = which(apply(v[b, , drop = F], 1, function(x) all(x == target)))
    if (length(ind) != 0)
    {
      rst = lapply(b[ind], function(x) c(hopeV, x))
    }
  }
  else
  {
    # cat("Find bounds input:\n")
    # cat(lb - 1L, ";", ub - 1L, "\n", sep = " ")
    # cat("Find bounds output:\n")


    bounds = findBoundsR(len = len, v = v, target = target,
                         me = 0L, lb = lb, ub = ub)


    # if(bounds$foundOrNot != 0L)
    #   cat(bounds$lb - 1L, ";", bounds$ub - 1L, "\n", sep = " ")
    # else cat("Contraction fails\n")


    if ( bounds$foundOrNot == 2L )
    {
      rst = list(c(hopeV, bounds$ub))
    }


    else if ( bounds$foundOrNot == 1L )
    {
      lb = bounds$lb
      ub = bounds$ub


      absDiff = ub - lb
      whichZero = which(absDiff == 0L)
      newHope = ub[whichZero]
      hopeV = c(hopeV, newHope)


      if ( length(whichZero) != 0 )
      {
        target = target - colSums(v[ub[whichZero], , drop = F])
        lb = lb[-whichZero]
        ub = ub[-whichZero]
        len = len - length(whichZero)
      }


      I = which.min(ub - lb)
      l = (lb[I] + ub[I]) %/% 2L
      u = l + 1L


      # The range at position I is [1, l]
      if(I != len) ubnew = c(pmin((l - I + 1L):l, ub[1:I]), ub[(I + 1L):len])
      else ubnew = pmin((l - I + 1L):l, ub[1:I])


      # The range at position I is [u, len]
      if(I != 1L) lbnew = c(lb[1:(I - 1L)], pmax(u:(u + len - I), lb[I:len]))
      else lbnew = pmax(u:(u + len - I), lb[I:len])


      # cat("I = ", I - 1L, ", len = ", len,
      #     ", Bounds after halving space:\n", sep = "")


      if ((I - 1L) <= (len - 1L) %/% 2L)
      {
        # cat(lb - 1L, ";", ubnew - 1L, "\n", sep = " ")
        # cat("\n")
        rstLeft = mflsssArbR(len, v, target, lb = lb, ub = ubnew, hopeV = hopeV)
        rstRight = mflsssArbR(len, v, target, lb = lbnew, ub = ub, hopeV = hopeV)

      }
      else
      {
        # cat(lbnew - 1L, ";", ub - 1L, "\n", sep = " ")
        # cat("\n")
        rstRight = mflsssArbR(len, v, target, lb = lbnew, ub = ub, hopeV = hopeV)
        rstLeft = mflsssArbR(len, v, target, lb = lb, ub = ubnew, hopeV = hopeV)
      }


      rst = c(rstLeft, rstRight)
    }


  }; rst
}


# Debug mflsssArbR
if(F)
{



  Rcpp::sourceCpp("src/mflsssArb.cpp", verbose = T)


  load("debug.Rdata")
  target = round(colSums(X[sol, , drop = F]), roundN)
  tmpX = round(X * 10 ^ roundN)
  tmpXmin = apply(tmpX, 2, function(x) min(x))
  tmpX = t(t(tmpX) - tmpXmin)
  tmpTarget = round(target * 10 ^ roundN)
  tmpTarget = round(tmpTarget - tmpXmin * len)
  tmpXrank = apply(tmpX, 2, function(x) rank(x, ties.method = "first"))
  colOdr = order(colSums(t(tmpXrank) %*% tmpXrank))
  tmpX = tmpX[, colOdr, drop = F]
  tmpTarget = as.integer(tmpTarget[colOdr])
  rowOdr = do.call(order, as.data.frame(tmpX[, ncol(tmpX):1, drop = F]))
  tmpX = matrix(as.integer(tmpX[rowOdr, , drop = F]), nrow = nrow(tmpX))
  tmpsol = sort(match(sol, rowOdr))


  rst = mflsssArbTest(len, v = tmpX, target = tmpTarget,
                      lb = 1:len, ub = (N - len + 1L):N, hopeV = integer(0))




}


Rcpp::sourceCpp("src/arbitraryDimFLSSStests/mflsssArb.cpp", verbose = T)
Rcpp::sourceCpp("src/arbitraryDimFLSSStests/findBound.cpp", verbose = T)
source("src/arbitraryDimFLSSStests/stringMatTo64bitIntMat.R")
source("R/numStrings.R")



# Tests.
if(F)
{

  # 2^64 = "18446744073709551616"
  #        "18446744073709551615"


  # Rcpp::sourceCpp("src/arbitraryDimFLSSStests/findBound.cpp", verbose = T)
  Rcpp::sourceCpp("src/mflsssArb.cpp", verbose = T)


  set.seed(42)


  options(scipen = 999)


  iter = 0L
  while(T){
    iter = iter + 1L
    # if (iter %% 10L == 0L) cat(iter, "")
    cat(iter, "")


    d = sample(5, 1) + 1L
    N = sample(30, 1)
    len = sample(N, 1)
    roundN = sample(5, 1) + 1L


    # d = 5L
    # len = 6L
    # N = 15L
    # roundN = 5L


    X = matrix(round(runif(d * N, -1, 1), roundN), nrow = N)
    sol = sort(sample(N, len))


    save(N, d, len, X, sol, roundN, file = "debug.Rdata")



    # load("debug.Rdata")
    target = round(colSums(X[sol, , drop = F]), roundN)
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


    # Rcpp::sourceCpp("src/arbitraryDimFLSSStests/mflsssArb.cpp", verbose = T)
    # sink("debug.txt")
    system.time({rstObjs = decomposeArbFLSSS(
      len = len, V = Xc, target = targetc, approxNinstance = sample(100, 1),
      maxCore = 10)})


    rst = unlist(lapply(rstObjs$arbFLSSSobjects, function(x)
      arbFLSSSobjRun(x, solutionNeed = 1e9, tlimit = 1e9)), recursive = F)
    rst = c(rstObjs$solutionsFound, rst)


    system.time({correctRst = arbFLSSS(
      len = len, V = Xc, target = targetc, sizeNeeded = 1e9,
      maxCore = 10, approxNinstance = 10, tlimit = 1e9)}) #; rst



    # correctRst = FLSSS::mFLSSSpar(
    #   maxCore = 10L, len = len, mV = tmpX, mTarget = tmpTarget,
    #   mME = numeric(ncol(tmpX)) + 0.3, solutionNeed = 1e9L,
    #   tlimit = 1e9, useBiSrchInFB = F, avgThreadLoad = 10L)


    # tmp = FLSSS::decomposeMflsss(
    #   len, tmpX, tmpTarget, numeric(ncol(tmpX)) + 0.3, solutionNeed = 1e9,
    #   approxNinstance = 5)


    rst = sort(unlist(lapply(rst, function(x) paste0(sort(x), collapse = "-"))))
    correctRst = sort(unlist(lapply(correctRst, function(x)
      paste0(sort(x), collapse = "-"))))


    if ( length(rst) != length(correctRst)) break
    if ( !all(unlist(rst) == unlist(correctRst))) break


  }






}




# Just a algo prototype for brute-force fixed-size subset sum.
if(F)
{


  reticulate::source_python('src/arbitraryDimFLSSStests/btsubsetsum.py')
  maxval = 300000L
  N = 10L
  len = 6L
  v = sort(sample(maxval, N))
  lb = (1:len) - 1L
  ub = (N - len):(N - 1L)
  rst = sort(btsubsetSum(v, lb, ub))
  tmp = sort(unique(sapply(1:100000, function(x) sum(sample(v, len)))))
  range(rst - tmp)




  iter = 0L
while(T){
  cat(iter, "")
  iter = iter + 1L

  N = 30L
  len = sample(N, 1)
  v = sort(sample(maxval, N))
  lb = sort(sample(N, len) - 1L)
  ub = sort(sample(N, len) - 1L)
  tmp = list(pmin(lb, ub), pmax(lb, ub))
  lb = tmp[[1]]; ub = tmp[[2]]
  rst = sort(btsubsetSum(v, lb, ub, uniqueOut = F))
  d = abs(NofSums(lb, ub) - length(rst)); d


  if (d != 0) break


}




}































