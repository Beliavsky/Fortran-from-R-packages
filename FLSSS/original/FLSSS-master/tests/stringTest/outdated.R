

# Test addPositiveStrings()
if (FALSE)
{


  x = sample(9, 20, replace = FALSE)
  x = sapply(x, function(x)
  {
    paste0(sample(0:9, x), collapse = '')
  })


  rst = addPositiveStrings(x)
  correctRst = sum(as.integer(x))
  as.numeric(rst) - correctRst


}


# Test subtract2positiveStrings()
if(FALSE)
{

  while (TRUE){
    x = paste0(sample(9, sample(10, 1), replace = TRUE) - 1L, collapse = '')
    y = paste0(sample(9, sample(10, 1), replace = TRUE) - 1L, collapse = '')
    z = subtract2positiveStrings(x, y)
    correctRst = as.numeric(x) - as.numeric(y)
    tmp = as.numeric(z) - correctRst
    if(abs(tmp) >= 0.9) break
  }


}



if(FALSE)
{


  NfracDigits = 7L
  while (TRUE) {

    s = round(runif(10, min = -1000, max = 1000), NfracDigits)


    rst = addNumStrings(as.character(s))


    correctRst = as.character(round(sum(as.numeric(s)), NfracDigits))



    if(rst != correctRst) break
  }


}


# Test test_intString2bitIntArray
if (FALSE)
{


  # 2^2048



  set.seed(123)


  # while(T){
  x = sample(2, 64, replace = TRUE) - 1L; x[1] = 0L
  y = sample(2, 64, replace = TRUE) - 1L; y[1] = 0L
  tmp = testMulHint(x, y)
  # if(tmp$rst != tmp$correctRst) break
  # }


  # testMulHintV(IntegerVector xv, IntegerVector yv)
  Rcpp::sourceCpp('src/intString2intArray.cpp', verbose = TRUE)
  set.seed(123)


  n = 2L
  while (TRUE){
    x = unlist(lapply(1:n, function(x) {
      x = sample(2, 64, replace = TRUE) - 1L; x[1] = 0L; x }))
    y = unlist(lapply(1:n, function(x) {
      x = sample(2, 64, replace = TRUE) - 1L; x[1] = 0L; x }))


    # tmp = testAddHintV(x, y)
    tmp = testSubHintV(x, y)


    if(tmp$rst != tmp$correctRst) break
  }


  if(FALSE)
  {


    load("tmp.Rdata")
    Rcpp::sourceCpp('src/intString2intArray.cpp')
    tmp = testMulHintV(x, y)


    source("R/numStrings.R")
    tmpf = function(i = 2048L)
    {
      rst = list("1")
      for(k in 1:(i - 1L))
      {
        x = rst[[length(rst)]]
        rst[[length(rst) + 1]] = addIntStrings(c(x, x))
      }
      unlist(rst)
    }
    power2string = tmpf()



  }



  n = 10L
  while (TRUE) {


    x = sample(2, 64 * n, replace = TRUE) - 1L
    y = sample(2, 64 * n, replace = TRUE) - 1L
    x[1] = 0L; y[1] = 0L


    # x = rep(1L, 64 * n)
    # y = rep(1L, 64 * n)
    # x[1] = 0L; y[1] = 0L

    # save(x, y, file = "tmp.Rdata")
    # tmp = testMulHintVbinInput(x, y)
    tmp = testAddHintVbinInput(x, y)
    tmp$rst == tmp$correctRst


    if(tmp$rst != tmp$correctRst) break
  }




  source("R/numStrings.R")
  tmpf = function(i = 2048L)
  {
    rst = list("1")
    for(k in 1:(i - 1L))
    {
      x = rst[[length(rst)]]
      rst[[length(rst) + 1]] = addIntStrings(c(x, x))
    }
    unlist(rst)
  }
  power2string = tmpf()


  tmpf = function(x)
  {
    a = which(rev(x) == 1L)
    if(length(a) == 0) return("0")
    addIntStrings(power2string[a])
  }
  set.seed(123)
  Rcpp::sourceCpp('src/intString2intArray.cpp', verbose = TRUE)
  while (TRUE){
    Ndigit = sample(300L, 1)
    s = paste0(c(sample(9, 1), sample(10, Ndigit - 1, replace = TRUE) - 1L), collapse = '')
    # s = '1234567890123'
    s2int2bin = testIntString2hugeInt32(s)
    s2int2bin2str = tmpf(s2int2bin)
    # s; s2int2bin2str
    if(s != s2int2bin2str) break
  }




  source("R/numStrings.R")
  tmpf = function(i = 4096L * 2L)
  {
    rst = list("1")
    for(k in 1:(i - 1L))
    {
      x = rst[[length(rst)]]
      rst[[length(rst) + 1]] = addIntStrings(c(x, x))
    }
    unlist(rst)
  }
  power2string = tmpf()


  binIntArray2intString = function(x)
  {
    a = which(rev(x) == 1L)
    if(length(a) == 0) return("0")
    addIntStrings(power2string[a])
  }


  mulIntString = function(x, y)
  {
    if(nchar(x) < nchar(y)) { tmp = x; x = y; y = tmp }
    y = as.integer(strsplit(y, '')[[1]])
    rst = list()
    for(i in length(y):1)
    {
      tmp = addIntStrings(rep(x, y[i]))
      rst[[length(rst) + 1]] =
        paste0(tmp, paste0(rep(0, length(y) - i), collapse = ''))
    }
    addIntStrings(unlist(rst))
  }


  set.seed(123)
  Rcpp::sourceCpp('src/intString2intArray.cpp', verbose = TRUE)


  while (TRUE){
    Ndigit = sample(1000L, 1)
    x = paste0(c(sample(9, 1), sample(10, Ndigit - 1, replace = TRUE) - 1L),
               collapse = '')
    Ndigit = sample(1000L, 1)
    y = paste0(c(sample(9, 1), sample(10, Ndigit - 1, replace = TRUE) - 1L),
               collapse = '')
    tmp = binIntArray2intString(testMulHintV32(x, y))


    tmpcorrect = mulIntString(x, y)
    tmp == tmpcorrect


    if(tmpcorrect != tmp) break
  }



}


# X is the input floating-point string matrix.
# len is the subset size
toIntStringMat = function(X, len, target)
{
  if(is.vector(X)) X = as.matrix(X)
  X = as.list(as.data.frame(X))
  target = as.list(target)
  largestSumVals = list()
  ranks = list()
  newtargets = list()
  for(i in 1:length(target))
  {
    tg = target[[i]]
    x = c(X[[i]], tg)
    validateNumStrings(x)
    x = asIntegerString(x)$s
    tg = x[length(x)]
    x = x[-length(x)]


    xNoSign = gsub('-', '', x)
    maxnchar = max(nchar(xNoSign))
    xNoSign = sapply(xNoSign, function(u)
    {
      zeros = paste0(rep('0', maxnchar - nchar(u)), collapse = '')
      paste0(zeros, u)
    })
    whichNeg = substr(x, 1, 1) == '-'
    xNoSign[whichNeg] = paste0('-', xNoSign[whichNeg])
    x = xNoSign


    id = order(x)
    ysorted = x[id]
    whichNeg = substr(ysorted, 1, 1) == '-'
    ysorted[whichNeg] = rev(ysorted[whichNeg])
    id[whichNeg] = rev(id[whichNeg])
    negYmin = ysorted[1]
    if(substr(negYmin, 1, 1) == '-')
      negYmin = substr(negYmin, 2, nchar(negYmin))
    else negYmin = paste0('-', negYmin)
    ysortedShifted = sapply(ysorted, function(u) addIntStrings(c(u, negYmin)))
    largestSubsetSum = addIntStrings(tail(ysortedShifted, n = len))
    addumForTarget = addIntStrings(rep(negYmin, len))
    tg = addIntStrings(c(tg, addumForTarget))
    if(substr(tg, 1, 1) == '-') return(list())
    newtargets[[i]] = tg
    largestSumVals[[i]] = largestSubsetSum
    X[[i]][id] = ysortedShifted
    rk = id
    rk[id] = 1:length(id)
    ranks[[i]] = rk
  }
  ranks = as.data.frame(ranks)
  names(ranks) = NULL
  rownames(ranks) = NULL
  list(X = as.data.frame(X, stringsAsFactors = FALSE), target = unlist(newtargets),
       ranks = ranks, largestSubsetSumVals = unlist(largestSumVals), len = len)
}


# Test toIntStringMat
if (FALSE)
{


  Rcpp::sourceCpp('src/validateStringInput.cpp', verbose = TRUE)


  d = 1L
  len = 15L
  size = 100L
  nd = 3L
  options(scipen = 999)


  while (TRUE){


  X = matrix(round(runif(d * size, -10, 10), nd), nrow = size)
  if(d == 1) target = round(sum(X[sample(size, len), ]), nd)
  else target = round(colSums(X[sample(size, len), ]), nd)
  Xc = matrix(as.character(X), nrow = size)
  targetc = as.character(target)
  rst = toIntStringMat(Xc, len, targetc)


  tmpcorrect = mapply(function(x, y)
  {
    y = as.integer(round((y - min(x) * len) * (10 ^ nd)))
    x = as.integer(round((x - min(x)) * (10 ^ nd)))
    rk = rank(x, ties.method = 'first')
    maxSubsetSum = sum(tail(sort(x), n = len))
    list(x = x, tg = y, rk = rk, maxSubsetSum = maxSubsetSum)
  }, as.data.frame(X), as.list(target), SIMPLIFY = FALSE)
  tmp = list()
  tmp$X = as.data.frame(lapply(tmpcorrect, function(x) x$x))
  tmp$target = unlist(lapply(tmpcorrect, function(x) x$tg))
  names(tmp$target) = NULL
  tmp$ranks = as.data.frame(lapply(tmpcorrect, function(x) x$rk))
  tmp$largestSubsetSumVals = unlist(lapply(tmpcorrect, function(x) x$maxSubsetSum))
  names(tmp$largestSubsetSumVals) = NULL
  tmp$len = len


  err = sum(abs(range(as.integer(unlist(rst$X)) - tmp$X))) +
    sum(abs(range(as.integer(rst$target) - tmp$target))) +
    sum(range(unlist(mapply(function(x, y) y[order(x)], rst$ranks, as.data.frame(X))) -
    unlist(mapply(function(x, y) y[order(x)], tmp$ranks, as.data.frame(X))))) +
    # sum(abs(range(unlist(rst$ranks) - tmp$ranks))) +
    sum(abs(range(unlist(as.integer(rst$largestSubsetSumVals)) -
                    tmp$largestSubsetSumVals))); err


  if(err != 0 ) break
}




}


# X is the return value from toIntStringMat
# list(X, target, ranks, largestSubsetSumVals, len)
# char dataframe, charvector, integer dataframe, charvector, integer
stringMatToBin = function(X, len, target)
{
  Y = toIntStringMat(X, len, target)
  corr = cor(Y$ranks)
  corrSum = colSums(corr)
  # odr = order(corrSum)
  # X = Y$X[, odr, drop = F]
  X = Y$X
  # target = Y$target[odr]
  target = Y$target
  # largestSubsetSumVals = Y$largestSubsetSumVals[odr]
  largestSubsetSumVals = Y$largestSubsetSumVals
  len = Y$len
  Xmaxcol = X[[which.max(corrSum)]]
  odr = order(nchar(Xmaxcol), Xmaxcol)
  # odr = order(X[[which.max(corrSum)]])
  X = X[odr, , drop = F]
  binaryValueIntList = list()
  targetValueIntList = list()
  largestSSintList = list()
  for(i in 1:length(X))
  {
    x = c(X[[i]], target[i], largestSubsetSumVals[i])
    binIntV = lapply(x, function(u) stringToBinIntVec(u))
    NlowZeros = min(unlist(lapply(binIntV, function(u) leftTrailingZeros(u))))
    if(NlowZeros != 0)
      binIntV = lapply(binIntV, function(u) u[-(1:NlowZeros)])
    Nslots = length(binIntV[[length(binIntV)]])
    binIntV = lapply(binIntV, function(u) c(u, integer(Nslots - length(u))))
    binaryValueIntList[[i]] = binIntV[1:(length(binIntV) - 2L)]
    targetValueIntList[[i]] = binIntV[[length(binIntV) - 1L]]
    largestSSintList[[i]] = binIntV[[length(binIntV)]]
  }
  Xbin = as.matrix(as.data.frame(lapply(1:length(binaryValueIntList[[1]]), function(i)
  {
    unlist(lapply(binaryValueIntList, function(u) u[[i]]))
  })))
  dimnames(Xbin) = NULL
  targetBin = unlist(targetValueIntList)
  largestSSbin = unlist(largestSSintList)
  names(targetBin) = NULL
  list(backOdr = odr, target = targetBin, X = Xbin, largestSSbin = largestSSbin)
}


# Test stringMatToBin
if(FALSE)
{


  Rcpp::sourceCpp('src/validateStringInput.cpp', verbose = TRUE)
  Rcpp::sourceCpp('src/intString2intArray.cpp', verbose = TRUE)
  source("R/numStrings.R")


  d = 2L
  len = 3L
  size = 5L
  nd = 2L
  options(scipen = 999)


  while (TRUE){


    X = matrix(round(runif(d * size, -10, 10), nd), nrow = size)
    if(d == 1) target = round(sum(X[sample(size, len), ]), nd) else
      target = round(colSums(X[sample(size, len), ]), nd)
    Xc = matrix(as.character(X), nrow = size)
    targetc = as.character(target)
    rst = toIntStringMat(Xc, len, targetc)
    tmp = stringMatToBin(Xc, len, targetc)


    tmp2 = AdjustStringVecTest(Xc, tsum = targetc, len)


    tmpcorrect = mapply(function(x, y)
    {
      y = as.integer(round((y - min(x) * len) * (10 ^ nd)))
      x = as.integer(round((x - min(x)) * (10 ^ nd)))
      rk = rank(x, ties.method = 'first')
      maxSubsetSum = sum(tail(sort(x), n = len))
      list(x = x, tg = y, rk = rk, maxSubsetSum = maxSubsetSum)
    }, as.data.frame(X), as.list(target), SIMPLIFY = FALSE)
    tmp = list()
    tmp$X = as.data.frame(lapply(tmpcorrect, function(x) x$x))
    tmp$target = unlist(lapply(tmpcorrect, function(x) x$tg))
    names(tmp$target) = NULL
    tmp$ranks = as.data.frame(lapply(tmpcorrect, function(x) x$rk))
    tmp$largestSubsetSumVals = unlist(lapply(tmpcorrect, function(x) x$maxSubsetSum))
    names(tmp$largestSubsetSumVals) = NULL
    tmp$len = len


    err = sum(abs(range(as.integer(unlist(rst$X)) - tmp$X))) +
      sum(abs(range(as.integer(rst$target) - tmp$target))) +
      sum(range(unlist(mapply(function(x, y) y[order(x)], rst$ranks, as.data.frame(X))) -
                  unlist(mapply(function(x, y) y[order(x)], tmp$ranks, as.data.frame(X))))) +
      # sum(abs(range(unlist(rst$ranks) - tmp$ranks))) +
      sum(abs(range(unlist(as.integer(rst$largestSubsetSumVals)) -
                      tmp$largestSubsetSumVals))); err


    if(err != 0 ) break
  }




}


