

# # Uncomment for Ubuntu
# setwd("/mnt/c/Users/i56087/Desktop/develop/develop036")


Rcpp::sourceCpp('src/arbitraryDimFLSSStests/intString2intArray.cpp', verbose = T)
source("R/numStrings.R")


binize = function(z)
{
  z = c(z, 0L)
  for(i in 1:(length(z) - 1L))
  {
    q = z[i] %/% 2L
    r = z[i] %% 2L
    z[i] = r
    z[i + 1] = z[i + 1] + q
  }
  z[1:max(which(z != 0))]
}


# Compute the correct result.
getCorrectRst = function(X, target, len, roundN)
{
  # orderCol = which.max(colSums(cor(apply(X, 2, function(x) rank(x, ties.method = "first")))))
  # odr = order(X[, orderCol])
  # X = as.data.frame(X[odr, ])


  orderCol = order(colSums(cor(apply(X, 2, function(x)
    rank(x, ties.method = "first")))))
  # odr = order(X[, orderCol])
  # X = as.data.frame(X[odr, ])
  X = as.data.frame(X[, orderCol])
  target = target[orderCol]


  target = as.list(target)
  for(i in 1:ncol(X))
  {
    minx = min(X[[i]])
    X[[i]] = X[[i]] - minx
    target[[i]] = target[[i]] - minx * len
  }


  scaler = 10 ^ roundN
  largestSS = list()
  Nbits = list()
  largestSSbit = list()
  for(i in 1:ncol(X))
  {
    X[[i]] = as.integer(round(X[[i]] * scaler))
    target[[i]] = as.integer(round(target[[i]] * scaler))
    if(target[[i]] < 0) return(list())
    largestSS[[i]] = sum(rev(sort(X[[i]]))[1:len])
    tmpS = as.integer(intToBits(largestSS[[i]]))
    Nbits[[i]] = tail(which(tmpS != 0L), n = 1)
    largestSSbit[[i]] = tmpS[1:Nbits[[i]]]
  }
  largestSSbit = unlist(largestSSbit)


  Xbit = list()
  targetBit = list()
  for(i in 1:length(X))
  {
    Xbit[[i]] = lapply(X[[i]], function(x)
    {
      a = integer(Nbits[[i]])
      b = as.integer(intToBits(x))
      l = min(length(a), length(b))
      a[1:l] = b[1:l]; a
    })


    a = integer(Nbits[[i]])
    b = as.integer(intToBits(target[[i]]))
    l = min(length(a), length(b))
    a[1:l] = b[1:l]
    targetBit[[i]] = a
  }
  targetBit = unlist(targetBit)


  Xbit = t(as.data.frame(lapply(1:length(Xbit[[1]]), function(i)
  {
    unlist(lapply(Xbit, function(x) x[[i]]))
  })))
  dimnames(Xbit) = NULL


  X = as.matrix(X); dimnames(X) = NULL
  target = unlist(target); dimnames(target) = NULL
  largestSS = unlist(largestSS); dimnames(largestSS) = NULL


  # list(XintBeforeOrder = X[order(odr), ], order = odr, Xint = X,
  #      targetInt = target, largestSSint = largestSS,
  #      Xbit = Xbit, targetBit = targetBit, largestSSbit = largestSSbit)


  odr = order(apply(Xbit, 1, function(x) paste0(rev(x), collapse = "")))
  XbitBeforeOrder = Xbit
  Xbit = Xbit[odr, ]


  list(targetInt = target, largestSSint = largestSS,
       XbitBeforeOrder = XbitBeforeOrder, targetBit = targetBit,
       largestSSbit = largestSSbit, Xbit = Xbit,
       order = odr, colOrder = orderCol)
}


intArrayToBits = function(x)
{
  y = as.integer(intToBits(x))
  whichNon0 = which(y != 0)
  if(length(whichNon0) == 0) return(0L)
  y[1:tail(whichNon0, n = 1)]
}


stringMatTo64bitIntMat = function(Xc, targetc, len, maxCore = 1)
{
  rst = stringMatTo64bitIntMatTest(Xc, targetc, len, maxCore)
  largestSSbit = intArrayToBits(rst$largestSSbit)
  l = length(largestSSbit)
  targetBit = intArrayToBits(rst$targetBit)
  targetBit = c(targetBit, integer(l - length(targetBit)))
  Xbit = t(as.data.frame(lapply(as.data.frame(t(rst$Xbit)), function(x)
  {
    y = intArrayToBits(x)
    c(y, integer(l - length(y)))
  })))
  dimnames(Xbit) = NULL
  list(order = rst$order, Xbit = Xbit, targetBit = targetBit,
       largestSSbit = largestSSbit, colOrder = rst$colOrder)
}




# Test body.
if(F)
{


  set.seed(123)


  k = 0L
  # while(T){
  k = k + 1L
  for (k in 1:10000){


    if(k %% 100L == 0L) cat(k, "")


    d = 10L
    len = 20L
    size = 100L
    roundN = 5L
    options(scipen = 999)
    X = matrix(round(runif(d * size, -10, 10), roundN), nrow = size)
    # save(X, roundN, size, len, d, file = "tmp2.Rdata")


    if(d == 1) target = round(sum(X[sample(size, len), ]), roundN) else
      target = round(colSums(X[sample(size, len), ]), roundN)
    Xc = matrix(as.character(X), nrow = size)
    targetc = as.character(target)


    correctRst = getCorrectRst(X, target, len, roundN)


    # rst = stringMatTo64bitIntMatTest(Xc, targetc, len, maxCore = 1)
    rst = stringMatTo64bitIntMat(Xc, targetc, len, maxCore = 10)


    err = sum(abs(unlist(correctRst[names(rst)]) - unlist(rst))) # err


    if(err != 0 & sum(abs(rst$order - correctRst$order)) == 0 &
       sum(abs(rst$colOrder - correctRst$colOrder)) == 0) break


  }


}































