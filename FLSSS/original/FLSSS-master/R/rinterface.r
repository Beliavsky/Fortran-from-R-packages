

premine <- function(len, v, target, ME)
{
  if(len == 0)
  {
    cat("Looping over subset size from 1 to the superset size is strongly recommended. Given zero subset size, setting a larger number of solutions in demand might be necessary to obtain the expected number of solutions.\n")
    return(NULL)
  }
  if(length(ME) == 1L) v = as.matrix(v)
  if(len == 1)
  {
    rst = as.list(which(apply(v, 1, function(x) all(abs(x - target) <= ME))))
    return(rst)
  }
  if(len == nrow(v) - 1L)
  {
    s = colSums(v)
    rst = as.list(which(apply(v, 1, function(x) all(abs(s - x - target) <= ME))))
    rst = lapply(rst, function(x) setdiff(1:nrow(v), x))
    return(rst)
  }
  if(len == nrow(v))
  {
    s = colSums(v)
    if(all(abs(s - target) <= ME)) return(list(1:nrow(v)))
    else return(list())
  }
  return(NULL)
}


# =============================================================================
# FLSSSvariableTree() will not be exposed and only serve research purpose.
# =============================================================================
FLSSSvariableTree <- function(len, v, target, ME, solutionNeed = 1L, LB = 1L : len, UB = (length(v) - len + 1L) : length(v), viaConjugate = FALSE, tlimit = 60, useBiSrchInFB = FALSE, useFloat = FALSE)
{
  premineRst = premine(len, v, target, ME)
  if(!is.null(premineRst)) return(premineRst)


  if(len == 0)
  {
    len = length(v)
    v = c(rep(0, len), v)
    vindex = c(rep(0L, len), 1L : len)
    sortOrder = order(v)
    v = v[sortOrder]
    vindex = vindex[sortOrder]
    rst = z_FLSSSvariableTree(len, v, target, ME, LB = 1L : len, UB = (length(v) - len + 1L) : length(v), solutionNeed, tlimit, useBiSrchInFB, useFloat)
    rst = unique(lapply(rst, function(x) sort(vindex[x][vindex[x] > 0L])))
    return(rst)
  }


  if(is.null(viaConjugate))
  {
    if(2L * len < length(v)) viaConjugate = T
  }


  if(viaConjugate)
  {
    target = sum(v) - target
    LBresv = LB
    LB = (1L : length(v))[-UB]
    UB = (1L : length(v))[-LBresv]
    len = length(v) - len
  }


  rst = z_FLSSSvariableTree(len, v, target, ME, LB, UB, solutionNeed, tlimit, useBiSrchInFB, useFloat)


  if(viaConjugate)
  {
    tmp = 1L : length(v)
    rst = lapply(rst, function(x) tmp[-x])
  }
  rst
}


# mV is the data matrix, each row is an observation
mFLSSSparVariableTree <- function(maxCore = 7L, len, mV, mTarget, mME, viaConjugate = NULL, solutionNeed = 1L, tlimit = 60, dl = ncol(mV), du = ncol(mV), randomizeTargetOrder = FALSE, useBiSrchInFB = FALSE, useFloat = FALSE)
{
  premineRst = premine(len, mV, mTarget, mME)
  if(!is.null(premineRst)) return(premineRst)


  # source("../src/legacy/rfuns.r")
  if(is.matrix(mV)) mV = as.data.frame(mV)


  dl = dl + 0L # materialize the values
  du = du + 0L


  fixedSize = T
  if(len == 0)
  {
    fixedSize = F
    len = nrow(mV)
    mV = as.data.frame(lapply(mV, function(x)c(rep(0, len), x)))
    randomizeTargetOrder = T
  }
  backout = function(mV, len){unique(lapply(mV, function(x) sort(x[x > len] - len)))}


  shouldConjugate = F
  if(is.null(viaConjugate) & 2L * len > nrow(mV)) shouldConjugate = T
  else if(!is.null(viaConjugate))
  {
    if(viaConjugate) shouldConjugate = T
  }


  if(shouldConjugate)
  {
    len = nrow(mV) - len
    mTarget = colSums(mV) - mTarget
  }


  if(ncol(mV) == 1L)
  {
    cat("Please call FLSSS for single dimensional set.\n")
    return(list())
  }


  LB = 1L : len
  UB = (nrow(mV) - len + 1L) : nrow(mV)


  # _______________________________________________________________________________________________
  # If mV are comonotonic
  # if(min(unlist(lapply(mV, function(x) min(diff(x))))) >= 0)
  # {
  #   d = ncol(mV)
  #   return(z_mFLSSScomo(len, mV, d, 0L, dl, d - du, du, mTarget, mME, LB, UB, solutionNeed, tlimit, useFloat, useBiSrchInFB))
  # }
  # _______________________________________________________________________________________________




  # _______________________________________________________________________________________________
  # z_viaLeadingCol <- function(len, mV, ME, target, randomOrderKeyTarget)
  info = z_viaLeadingCol(len, mV, mME, mTarget, randomizeTargetOrder, dl, du, F)
  if(is.character(info)) return(info)


  # _______________________________________________________________________________________________
  # Find the leading column that most correlates the rest columns, order other columns by the leading column. In the recent research, this substantially accelerates the speed.
  leadingCol = which.max(colSums(cor(mV, method = "spearman")))
  leadingColOrder = order(mV[[leadingCol]])
  mV = mV[leadingColOrder, ]
  info2 = z_adhereIndexCol(len, mV, mME, mTarget, randomizeTargetOrder, dl, du, F)
  if(!is.list(info2)) return(info2)


  # info = info2
  # It is mythical for now that, even length(info$keyTarget) > length(info2$keyTarget), adding index as the key column still makes things faster... Interesting, so I decided to do "/2", to encourage adding the index column.
  if(is.null(info) | length(info$keyTarget) > length(info2$keyTarget) / 2) info = info2
  # if(verbose) cat("Final atom tasks = ", length(info$keyTarget), "\n")


  # tlimit = tlimit * maxCore
  d = ncol(info$v)
  dl = info$dldu[1]
  du = info$dldu[2]
  zeroBasedKeyInd = info$zeroBasedKeyInd
  rst = z_mFLSSSvariableTree(maxCore, len, info$v, d, 0, dl, d - du, du, zeroBasedKeyInd, info$originalTarget, info$keyTarget, info$scaleFactor, info$ME, LB, UB, solutionNeed, tlimit, useFloat, useBiSrchInFB)
  rst = lapply(rst, function(x) leadingColOrder[x])
  if(shouldConjugate)
  {
    tmp = 1L : nrow(mV)
    rst = lapply(rst, function(x) tmp[-x])
  }
  if(!fixedSize) rst = backout(rst, len)
  rst
  # list(solution = rst, memoryImage = NA)
}


# =============================================================================
# mFLSSSparVariableTree() will not be exposed and only serve research purpose.
# =============================================================================




FLSSS <- function(len, v, target, ME, solutionNeed = 1L, LB = 1:len, UB = (length(v) - len + 1L):length(v), viaConjugate = FALSE, tlimit = 60, useBiSrchInFB = FALSE, NfractionDigits = Inf)
{
  valtype = "int"
  if(is.infinite(NfractionDigits)) valtype = "double"
  else
  {
    scaler = 10 ^ NfractionDigits
    v = as.integer(round(v * scaler))
    target = as.integer(round(target * scaler))
    ME = as.integer(round(ME * scaler))
  }


  premineRst = premine(len, v, target, ME)
  if(!is.null(premineRst)) return(premineRst)


  if(len == 0)
  {
    len = length(v)
    v = c(rep(0, len), v)
    vindex = c(rep(0L, len), 1L : len)
    sortOrder = order(v)
    v = v[sortOrder]
    vindex = vindex[sortOrder]
    if(valtype == "double")
    {
      target = target / ME
      v = v / ME
      ME = 1
    }
    rst = z_FLSSS(len, v, target, ME, LB = 1:len, UB = (length(v) - len + 1L):length(v), solutionNeed, tlimit, useBiSrchInFB, valtype)
    rst = unique(lapply(rst, function(x) sort(vindex[x][vindex[x] > 0L])))
    return(rst)
  }


  if(valtype == "double")
  {
    target = target / ME
    v = v / ME
    ME = 1
  }


  if(is.null(viaConjugate))
  {
    if(2L * len < length(v)) viaConjugate = T
  }


  if(viaConjugate)
  {
    target = sum(v) - target
    LBresv = LB
    LB = (1:length(v))[-UB]
    UB = (1:length(v))[-LBresv]
    len = length(v) - len
  }


  rst = z_FLSSS(len, v, target, ME, LB, UB, solutionNeed, tlimit, useBiSrchInFB, valtype)


  if(viaConjugate)
  {
    tmp = 1:length(v)
    rst = lapply(rst, function(x) tmp[-x])
  }
  rst
}








# -------------------------------------------------------------------------------------------------
FLSSSmultiset <- function(len, buckets, target, ME, solutionNeed = 1L, tlimit = 60, useBiSrchInFB = FALSE, NfractionDigits = Inf)
{


  if(length(buckets) == 1L) return("buckets size equals 1. Please call FLSSS().")
  if(min(len) <= 0L) return("len must be positive.")


  bucketsOrder = order(unlist(lapply(buckets, function(x) max(x))))
  buckets = buckets[bucketsOrder]
  len = len[bucketsOrder]
  orderInBucket = lapply(buckets, function(x) order(x))


  buckets = mapply(function(x, y) x[y], buckets, orderInBucket, SIMPLIFY = F)


  bucketsSize = unlist(lapply(buckets, function(x) length(x)))
  ub = cumsum(bucketsSize)
  lb = ub - bucketsSize + 1L


  LB = integer(sum(len))
  UB = LB
  whichBucket = UB
  i = 1L
  j = 1L
  for(j in 1L : length(len))
  {
    tmp = i : (i + len[j] - 1L)
    LB[tmp] = lb[j] : (lb[j] + len[j] - 1L)
    UB[tmp] = (ub[j] - len[j] + 1L) : ub[j]
    whichBucket[tmp] = j
    i = i + len[j]
  }


  for(i in 2:length(buckets))
  {
    target = target + (buckets[[i - 1L]][length(buckets[[i - 1L]])] - buckets[[i]][1]) * len[i]
    buckets[[i]] = (buckets[[i - 1L]][length(buckets[[i - 1L]])] - buckets[[i]][1]) + buckets[[i]]
  }
  v = unlist(buckets)


  valtype = "int"
  if(is.infinite(NfractionDigits)) valtype = "double"
  else
  {
    scaler = 10 ^ NfractionDigits
    v = as.integer(round(v * scaler))
    target = as.integer(round(target * scaler))
    ME = as.integer(round(ME * scaler))
  }


  result = z_FLSSS(length(LB), v, target, ME, LB, UB, solutionNeed, tlimit, useBiSrchInFB, valtype)


  result = lapply(result, function(x)
  {
    tmp = stats::aggregate(list(sort(x)), list(whichBucket), function(t) t, simplify = F)[[2]]
    tmp = mapply(function(o, a, b) o[a - b + 1L], orderInBucket, tmp, lb, SIMPLIFY = F)
    rst = vector("list", length(tmp))
    rst[bucketsOrder] = tmp
    rst
  })


  result
}








# -------------------------------------------------------------------------------------------------
# mV is the data matrix, each row is an observation
mFLSSSparImposeBounds <- function(maxCore = 7L, len, mV, mTarget, mME, LB = 1L : len, UB = (nrow(mV) - len + 1L) : nrow(mV), solutionNeed = 1L, tlimit = 60, dl = ncol(mV), du = ncol(mV), targetsOrder = NULL, useBiSrchInFB = FALSE, avgThreadLoad = 8L)
{
  premineRst = premine(len, mV, mTarget, mME)
  if(!is.null(premineRst)) return(premineRst)


  if(is.data.frame(mV)) mV = as.matrix(mV)


  d = ncol(mV)
  dlst = 0L
  dl = dl + 0L # materialize the values
  dust = d - du
  du = du + 0L
  if(d == 1)
  {
    cat("Please call FLSSS for single dimensional set.\n")
    return(list())
  }


  fixedSize = T
  if(len == 0)
  {
    fixedSize = F
    len = nrow(mV)
    mV = rbind(matrix(numeric(len * d), ncol = d), mV)
  }
  backout = function(rst, len){unique(lapply(rst, function(x) sort(x[x > len] - len)))}


  shouldConjugate = F
  if(2L * len > nrow(mV)) shouldConjugate = T


  if(shouldConjugate)
  {
    len = nrow(mV) - len
    mTarget = colSums(mV) - mTarget
    tmp = dlst; dlst = dust; dust = tmp
    tmp = dl; dl = du; du = tmp
  }


  # _______________________________________________________________________________________________
  # If mV are comonotonic
  if(min(unlist(apply(mV, 2, function(x) min(diff(x))))) >= 0)
  {
    d = ncol(mV)
    return(z_mFLSSScomoPar(maxCore, len, mV, numeric(0), d, dlst, dl, dust, du, mTarget, mME, LB, UB, solutionNeed, tlimit, useBiSrchInFB, avgThreadLoad))
  }
  # _______________________________________________________________________________________________


  # generate target matrix, each column is a d-dimensional target vector
  tmp = z_mTargetMatNoKey(len, mV, mTarget, mME, targetsOrder, dl, du)
  mV = tmp$mV
  targetMat = tmp$targetMat
  mME = tmp$mME
  d = tmp$d
  dl = tmp$dl
  du = tmp$du
  leadingCol = tmp$keyColInd
  rm(tmp); gc()


  # tlimit = tlimit * maxCore
  # print(len); print(mV); print(d); print(dlst); print(dust); print(du); # print(targetMat);
  # print(mME)
  rst = z_mFLSSS(maxCore, len, mV, numeric(0), d, dlst, dl, dust, du, targetMat, mME, LB, UB, solutionNeed, tlimit, useBiSrchInFB)


  if(shouldConjugate)
  {
    tmp = 1L : nrow(mV)
    rst = lapply(rst, function(x) tmp[-x])
  }
  if(!fixedSize) rst = backout(rst, len)
  rst
}







# target
# mV is the data matrix, each row is an observation
mFLSSSpar <- function(maxCore = 7L, len, mV, mTarget, mME, solutionNeed = 1L, tlimit = 60, dl = ncol(mV), du = ncol(mV), useBiSrchInFB = FALSE, avgThreadLoad = 8L)
{
  premineRst = premine(len, mV, mTarget, mME)
  if(!is.null(premineRst)) return(premineRst)


  if(is.data.frame(mV)) mV = as.matrix(mV)


  d = ncol(mV)
  dlst = 0L
  dl = dl + 0L # materialize the values
  dust = d - du
  du = du + 0L
  if(d == 1)
  {
    cat("Please call FLSSS for single dimensional set.\n")
    return(list())
  }


  # Make them all nonnegative, but only for fixed size subset sum.
  if(len != 0)
  {
    minV = apply(mV, 2L, function(x) min(x))
    mV = apply(mV, 2L, function(x) x - min(x))
    mTarget = mTarget - len * minV
  }


  fixedSize = T
  if(len == 0)
  {
    fixedSize = F
    len = nrow(mV)
    mV = rbind(matrix(0, nrow = len, ncol = d), mV)
  }
  backout = function(rst, len) { unique(lapply(rst, function(x) sort(x[x > len] - len)) ) }


  shouldConjugate = F
  if(2L * len > nrow(mV)) shouldConjugate = T


  if(shouldConjugate)
  {
    len = nrow(mV) - len
    mTarget = colSums(mV) - mTarget
    tmp = dlst; dlst = dust; dust = tmp
    tmp = dl; dl = du; du = tmp
  }


  LB = 1L : len
  UB = (nrow(mV) - len + 1L) : nrow(mV)


  # Find the leading column that most correlates the rest columns, order other columns by the leading column. Recent simulations show it substantially accelerates the speed.
  corMat = cor(mV, method = "spearman")
  corMat[is.na(corMat)] = 0
  leadingCol = which.max(colSums(corMat))
  leadingColOrder = order(mV[, leadingCol])
  mV = mV[leadingColOrder, ]


  # _______________________________________________________________________________________________
  # If mV are comonotonic
  if(min(unlist(apply(mV, 2, function(x) min(diff(x))))) >= 0)
  {
    d = ncol(mV)
    rst = z_mFLSSScomoPar(maxCore, len, mV, numeric(0), d, dlst, dl, dust, du, mTarget, mME, LB, UB, solutionNeed, tlimit, useBiSrchInFB, avgThreadLoad)
    rst = lapply(rst, function(x) leadingColOrder[x])
    if(shouldConjugate)
    {
      tmp = 1L : nrow(mV)
      rst = lapply(rst, function(x) tmp[-x])
    }
    if(!fixedSize) rst = backout(rst, len)
    return(rst)
  }
  # _______________________________________________________________________________________________


  # generate target matrix, each column is a d-dimensional target vector
  tmp = z_mTargetMat(len, mV, mTarget, mME, leadingCol, dl, du)
  mV = tmp$mV
  targetMat = tmp$targetMat
  mME = tmp$mME
  d = tmp$d
  dl = tmp$dl
  du = tmp$du
  leadingCol = tmp$keyColInd
  rm(tmp); gc()


  dimnames(mV) = NULL


  # parameterList = list(maxCore = maxCore, len = len, mV = mV, d = d, targetMat = targetMat, mME = mME, LB = LB, UB = UB, solutionNeed = solutionNeed, tlimit = tlimit)
  # save(parameterList, file = "parameterList.Rdata")


  rst = z_mFLSSS(maxCore, len, mV, numeric(0), d, dlst, dl, dust, du, targetMat, mME, LB, UB, solutionNeed, tlimit, useBiSrchInFB)


  rst = lapply(rst, function(x) leadingColOrder[x])
  # print(str(rst))
  if(shouldConjugate)
  {
    tmp = 1L : nrow(mV)
    rst = lapply(rst, function(x) tmp[-x])
  }


  if(!fixedSize) rst = backout(rst, len)
  rst
}




# target
# mV is the data matrix, each row is an observation
decomposeMflsss <- function(len, mV, mTarget, mME, solutionNeed = 1L, dl = ncol(mV), du = ncol(mV), useBiSrchInFB = FALSE, approxNinstance = 50000L)
{
  premineRst = premine(len, mV, mTarget, mME)
  if(!is.null(premineRst)) return(list(mflsssObjects = list(), solutionsFound = premineRst))


  if(is.data.frame(mV)) mV = as.matrix(mV)


  d = ncol(mV)
  dlst = 0L
  dl = dl + 0L # materialize the values
  dust = d - du
  du = du + 0L
  if(d == 1)
  {
    cat("Please call FLSSS for single dimensional set.\n")
    return(list())
  }


  # Make them all nonnegative
  if(len != 0)
  {
    minV = apply(mV, 2L, function(x) min(x))
    mV = apply(mV, 2L, function(x) x - min(x))
    mTarget = mTarget - len * minV
  }


  fixedSize = T
  if(len == 0)
  {
    fixedSize = F
    len = nrow(mV)
    mV = rbind(matrix(numeric(len * d), ncol = d), mV)
  }
  backout = function(rst, len) { unique(lapply(rst, function(x) sort(x[x > len] - len))) }


  shouldConjugate = F
  if(2L * len > nrow(mV)) shouldConjugate = T


  if(shouldConjugate)
  {
    len = nrow(mV) - len
    mTarget = colSums(mV) - mTarget
    tmp = dlst; dlst = dust; dust = tmp
    tmp = dl; dl = du; du = tmp
  }


  LB = 1L : len
  UB = (nrow(mV) - len + 1L) : nrow(mV)


  # Find the leading column that most correlates the rest columns, order other columns by the leading column. Recent simulations show it substantially accelerates the speed.
  corMat = cor(mV, method = "spearman")
  corMat[is.na(corMat)] = 0
  leadingCol = which.max(colSums(corMat))
  leadingColOrder = order(mV[, leadingCol])
  mV = mV[leadingColOrder, ]


  maskV = numeric(0)
  produceImage = function()
  {
    dimnames(mV) = NULL
    rst = z_mFLSSSimage(len, mV, maskV, d, dlst, dl, dust, du, as.matrix(mTarget), mME, LB, UB, solutionNeed, useBiSrchInFB, approxNinstance)
    rst$mflsssObjects = lapply(rst$mflsssObjects, function(x)
    {
      c(x, list(leadingColOrder = leadingColOrder), list(shouldConjugate = shouldConjugate), list(fixedSize = fixedSize), list(maskV = maskV), list(len = len), list(backout = backout))
    })
    if(length(rst$solutionsFound) != 0)
    {
      tmprst = rst$solutionsFound
      tmprst = lapply(tmprst, function(x) leadingColOrder[x])
      if(shouldConjugate)
      {
        tmp = 1L : nrow(mV)
        tmprst = lapply(tmprst, function(x) tmp[-x])
      }
      if(!fixedSize) tmprst = backout(tmprst, len)
      rst$solutionsFound = tmprst
    }; rst
  }


  # _______________________________________________________________________________________________
  # If mV are comonotonic
  if(min(unlist(apply(mV, 2, function(x) min(diff(x))))) >= 0)
  {
    d = ncol(mV)
    rst = produceImage()
    return(rst)
  }
  # _______________________________________________________________________________________________


  # Generate target matrix, each column is a d-dimensional target vector
  tmp = z_mTargetMat(len, mV, mTarget, mME, leadingCol, dl, du)
  mV = tmp$mV
  mTarget = tmp$targetMat
  mME = tmp$mME
  d = tmp$d
  dl = tmp$dl
  du = tmp$du
  leadingCol = tmp$keyColInd
  rm(tmp); gc()


  dimnames(mV) = NULL
  rst = produceImage(); rst
}




mFLSSSobjRun <- function(mflsssObj, solutionNeed = 1, tlimit = 60)
{
  rst = z_mFLSSSimport(mflsssObj, solutionNeed, tlimit)
  leadingColOrder = mflsssObj$leadingColOrder
  shouldConjugate = mflsssObj$shouldConjugate
  fixedSize = mflsssObj$fixedSize
  rst = lapply(rst, function(x) leadingColOrder[x])
  mV = mflsssObj$vr
  len = mflsssObj$len
  backout = mflsssObj$backout
  if(shouldConjugate)
  {
    tmp = 1:nrow(mV)
    rst = lapply(rst, function(x) tmp[-x])
  }
  if(!fixedSize) rst = backout(rst, len)
  rst
}




# -------------------------------------------------------------------------------------------------
# mV is the data matrix, each row is an observation
mFLSSSparImposeBoundsIntegerized <- function(maxCore = 7L, len, mV, mTarget, mME, LB = 1L : len, UB = (nrow(mV) - len + 1L) : nrow(mV), solutionNeed = 1L, precisionLevel = integer(ncol(mV)), returnBeforeMining = FALSE, tlimit = 60, dl = ncol(mV), du = ncol(mV), targetsOrder = NULL, useBiSrchInFB = FALSE, avgThreadLoad = 8L, verbose = TRUE)
{
  premineRst = premine(len, mV, mTarget, mME)
  if(!is.null(premineRst)) return(premineRst)


  if(.Machine$sizeof.pointer == 4L)
  {
    message("32-bit architecture unsupported")
    return()
  }


  if(is.data.frame(mV)) mV = as.matrix(mV)
  N = nrow(mV)
  d = ncol(mV)
  dlst = 0L
  dl = dl + 0L # materialize values
  dust = d - du
  du = du + 0L


  if(d == 1L)
  {
    cat("Please call FLSSS for single dimensional set.\n")
    return(list())
  }


  # Make them all nonnegative
  if(len != 0)
  {
    minV = apply(mV, 2L, function(x) min(x))
    mV = apply(mV, 2L, function(x) x - min(x))
    mTarget = mTarget - len * minV
  }


  # up-bound the subset sum range, just in case
  {
    lb = mTarget - mME
    ub = mTarget + mME
    leastSsum = apply(mV, 2L, function(x) sum(sort(x)[1L : len]))
    mostSsum = apply(mV, 2L, function(x) sum(sort(x)[(N - len + 1L) : N]))
    tmpGap = (mostSsum - leastSsum) * (1 / 64)
    lb = pmax(leastSsum - tmpGap, lb, numeric(length(lb)))
    ub = pmin(mostSsum + tmpGap, ub)
    mTarget = (lb + ub) / 2L
    mME = (ub - lb) / 2L
  }


  fixedSize = T
  if(len == 0L)
  {
    fixedSize = F
    len = nrow(mV)
    mV = rbind(matrix(numeric(len * d), ncol = d), mV)
  }
  backout = function(rst, len){unique(lapply(rst, function(x) sort(x[x > len] - len)))}


  shouldConjugate = F
  if(2L * len > nrow(mV)) shouldConjugate = T


  if(shouldConjugate)
  {
    len = nrow(mV) - len
    mTarget = colSums(mV) - mTarget
    tmp = dlst; dlst = dust; dust = tmp
    tmp = dl; dl = du; du = tmp
    LB = sort(setdiff(1L : nrow(mV), LB))
    UB = sort(setdiff(1L : nrow(mV), UB))
  }


  # integerize
  {
    tmp = z_integerize(len, mV, mTarget, mME, precisionLevel)
    mV = tmp$integerized
    mTarget = tmp$target
    mME = tmp$ME
  }


  {
    dimnames(mV) = NULL; dimnames(mTarget) = NULL; dimnames(mME) = NULL
    INT = list(mV = mV, mTarget = mTarget, mME = mME)
  }


  # if(returnBeforeMining) return(INT)


  # _______________________________________________________________________________________________
  # If mV are comonotonic
  if(min(unlist(apply(mV, 2L, function(x) min(diff(x))))) >= 0L)
  {
    tmp = z_filterTargetFindLargestMagnitude(len, mV, as.matrix(mTarget), mME)
    tmp = z_crunchIntegers(len, mV, tmp$targetMat, mME, dlst, dl, dust, du, tmp$maxMag)
    mV = tmp$mV
    mTarget = as.numeric(tmp$targetMat)
    mME = tmp$mME
    maskV = tmp$maskV


    compressedDim = ncol(mV)
    # if(returnBeforeMining) return(c(INT, list(compressedDim = compressedDim)))
    if(verbose) cat("Dimensionality reduced from", ncol(INT$mV), "to", ncol(mV), "\n")
    if(returnBeforeMining) return(list(solution = list(), INT = c(INT, list(compressedDim = compressedDim))))


    rst = z_mFLSSScomoPar(maxCore, len, mV, maskV, d, dlst, dl, dust, du, mTarget, mME, LB, UB, solutionNeed, tlimit, useBiSrchInFB, avgThreadLoad)
    if(shouldConjugate)
    {
      tmp = 1L : nrow(mV)
      rst = lapply(rst, function(x) tmp[-x])
    }
    if(!fixedSize) rst = backout(rst, len)
    # return(list(solution = rst, INT = INT))
    return(list(solution = rst, INT = c(INT, list(compressedDim = compressedDim))))
  }
  # _______________________________________________________________________________________________


  # generate target matrix, each column is a d-dimensional target vector
  {
    tmp = z_mTargetMatNoKeyINT(len, mV, mTarget, mME, targetsOrder, dl, du)
    mV = tmp$mV
    targetMat = tmp$targetMat
    mME = tmp$mME
    d = tmp$d
    dl = tmp$dl
    du = tmp$du
    leadingCol = tmp$keyColInd
    rm(tmp); gc()
  }


  dimnames(mV) = NULL


  # subtract minima again
  {
    minV = mV[1, ]
    mV = apply(mV, 2, function(x) x - x[1])
    tmp = minV * as.integer(len)
    targetMat = apply(targetMat, 2, function(x) x - tmp)
    dimnames(mV) = NULL
    dimnames(targetMat) = NULL
    dimnames(mME) = NULL
  }


  targetMatAndMaxMag = z_filterTargetFindLargestMagnitude(len, mV, targetMat, mME)
  targetMat = targetMatAndMaxMag$targetMat
  maxMag = targetMatAndMaxMag$maxMag


  # crunch integers
  {
    tmp = z_crunchIntegers(len, mV, targetMat, mME, dlst, dl, dust, du, maxMag)
    mV = tmp$mV
    targetMat = tmp$targetMat
    mME = tmp$mME
    maskV = tmp$maskV
    d = length(maskV)
    dlst = tmp$dim[1]; dl = tmp$dim[2]; dust = tmp$dim[3]; du = tmp$dim[4]
  }


  compressedDim = ncol(mV)
  if(returnBeforeMining) # return(c(INT, list(compressedDim = ncol(mV))))
    return(list(solution = list(), INT = c(INT, list(compressedDim = compressedDim))))
  if(verbose) cat("Dimensionality reduced from", ncol(INT$mV) + 1L, "to", ncol(mV), "\n")


  # tlimit = tlimit * maxCore
  rst = z_mFLSSS(maxCore, len, mV, maskV, d, dlst, dl, dust, du, targetMat, mME, LB, UB, solutionNeed, tlimit, useBiSrchInFB)


  if(shouldConjugate)
  {
    tmp = 1L : nrow(mV)
    rst = lapply(rst, function(x) tmp[-x])
  }
  if(!fixedSize) rst = backout(rst, len)
  # list(solution = rst, INT = INT)
  list(solution = rst, INT = c(INT, list(compressedDim = compressedDim)))
}




# target
# mV is the data matrix, each row is an observation
mFLSSSparIntegerized <- function(maxCore = 7L, len, mV, mTarget, mME, solutionNeed = 1L, precisionLevel = integer(ncol(mV)), returnBeforeMining = FALSE, tlimit = 60, dl = ncol(mV), du = ncol(mV), useBiSrchInFB = FALSE, avgThreadLoad = 8L, verbose = TRUE)
{
  premineRst = premine(len, mV, mTarget, mME)
  if(!is.null(premineRst)) return(premineRst)


  if(.Machine$sizeof.pointer == 4L)
  {
    message("32-bit architecture unsupported")
    return()
  }


  if(is.data.frame(mV)) mV = as.matrix(mV)
  len = as.integer(len)
  N = nrow(mV)
  d = ncol(mV)
  dlst = 0L
  dl = dl + 0L # materialize values
  dust = d - du
  du = du + 0L


  if(d == 1)
  {
    cat("Please call FLSSS for single dimensional set.\n")
    return(list())
  }


  # Make them all nonnegative
  if(len != 0)
  {
    minV = apply(mV, 2L, function(x) min(x))
    mV = apply(mV, 2L, function(x) x - min(x))
    mTarget = mTarget - len * minV
  }


  # up-bound the subset sum range, just in case
  {
    lb = mTarget - mME
    ub = mTarget + mME
    leastSsum = apply(mV, 2L, function(x) sum(sort(x)[1L : len]))
    mostSsum = apply(mV, 2L, function(x) sum(sort(x)[(N - len + 1L) : N]))
    tmpGap = mostSsum - leastSsum
    lb = pmax(leastSsum - tmpGap * (1 / 8), lb, numeric(length(lb)))
    ub = pmin(mostSsum + tmpGap * (1 / 8), ub)
    mTarget = (lb + ub) / 2
    mME = (ub - lb) / 2
  }


  fixedSize = T
  if(len == 0)
  {
    fixedSize = F
    len = nrow(mV)
    mV = rbind(matrix(numeric(len * d), ncol = d), mV)
  }
  backout = function(rst, len){unique(lapply(rst, function(x) sort(x[x > len] - len)))}


  shouldConjugate = F
  if(2L * len > nrow(mV)) shouldConjugate = T


  if(shouldConjugate)
  {
    len = nrow(mV) - len
    mTarget = colSums(mV) - mTarget
    tmp = dlst; dlst = dust; dust = tmp
    tmp = dl; dl = du; du = tmp
  }


  # integerize
  {
    tmp = z_integerize(len, mV, mTarget, mME, precisionLevel)
    mV = tmp$integerized
    minV = apply(mV, 2, function(x) min(x))
    mV = apply(mV, 2, function(x) x - min(x)) # zeroing minima again
    mTarget = tmp$target - len * minV
    mME = tmp$ME
  }


  {
    dimnames(mV) = NULL; dimnames(mTarget) = NULL; dimnames(mME) = NULL
    INT = list(mV = mV, mTarget = mTarget, mME = mME)
  }


  LB = 1L : len
  UB = (nrow(mV) - len + 1L) : nrow(mV)


  # Find the leading column that most correlates the rest columns, order other columns by the leading column. Recent simulations show it substantially accelerates the speed.
  corMat = cor(mV, method = "spearman")
  corMat[is.na(corMat)] = 0
  leadingCol = which.max(colSums(corMat))
  leadingColOrder = order(mV[, leadingCol])
  mV = mV[leadingColOrder, ]


  # _______________________________________________________________________________________________
  # If mV are comonotonic
  if(min(unlist(apply(mV, 2, function(x) min(diff(x))))) >= 0)
  {
    tmp = z_filterTargetFindLargestMagnitude(len, mV, as.matrix(mTarget), mME)
    tmp = z_crunchIntegers(len, mV, tmp$targetMat, mME, dlst, dl, dust, du, tmp$maxMag)
    mV = tmp$mV
    mTarget = as.numeric(tmp$targetMat)
    mME = tmp$mME
    maskV = tmp$maskV
    d = length(maskV)
    dlst = tmp$dim[1]; dl = tmp$dim[2]; dust = tmp$dim[3]; du = tmp$dim[4]


    if(verbose) cat("Dimensionality reduced from", ncol(INT$mV), "to", ncol(mV), "\n")
    if(returnBeforeMining) # return(c(INT, list(compressedDim = ncol(mV))))
      return(list(solution = list(), INT = c(INT, list(compressedDim = compressedDim))))


    rst = z_mFLSSScomoPar(maxCore, len, mV, maskV, d, dlst, dl, dust, du, mTarget, mME, LB, UB, solutionNeed, tlimit, useBiSrchInFB, avgThreadLoad)
    rst = lapply(rst, function(x) leadingColOrder[x])
    if(shouldConjugate)
    {
      tmp = 1L : nrow(mV)
      rst = lapply(rst, function(x) tmp[-x])
    }
    if(!fixedSize) rst = backout(rst, len)
    # return(solution = rst, INT = INT)
    return(list(solution = rst, INT = c(INT, list(compressedDim = compressedDim))))
  }
  # _______________________________________________________________________________________________


  # generate target matrix, each column is a d-dimensional target vector
  {
    tmp = z_mTargetMatINT(len, mV, mTarget, mME, leadingCol, dl, du)
    mV = tmp$mV
    targetMat = tmp$targetMat
    mME = tmp$mME
    d = tmp$d
    dl = tmp$dl
    du = tmp$du
    leadingCol = tmp$keyColInd
    rm(tmp); gc()
  }


  # subtract minima again
  {
    minV = mV[1, ]
    mV = apply(mV, 2, function(x) x - x[1])
    targetMat = targetMat - minV * as.integer(len)
    dimnames(mV) = NULL
    dimnames(targetMat) = NULL
    dimnames(mME) = NULL
  }


  targetMatAndMaxMag = z_filterTargetFindLargestMagnitude(len, mV, targetMat, mME)
  targetMat = targetMatAndMaxMag$targetMat
  maxMag = targetMatAndMaxMag$maxMag


  # crunch integers
  {
    tmp = z_crunchIntegers(len, mV, targetMat, mME, dlst, dl, dust, du, maxMag)
    mV = tmp$mV
    targetMat = tmp$targetMat
    mME = tmp$mME
    maskV = tmp$maskV
    d = length(maskV)
    dlst = tmp$dim[1]; dl = tmp$dim[2]; dust = tmp$dim[3]; du = tmp$dim[4]
  }


  compressedDim = ncol(mV)
  if(returnBeforeMining)
    return(list(solution = list(), INT = c(INT, list(compressedDim = compressedDim))))
  if(verbose) cat("Dimensionality reduced from", ncol(INT$mV) + 1L, "to", ncol(mV), "\n")


  # tlimit = tlimit * maxCore
  rst = z_mFLSSS(maxCore, len, mV, maskV, d, dlst, dl, dust, du, targetMat, mME, LB, UB, solutionNeed, tlimit, useBiSrchInFB)


  rst = lapply(rst, function(x) leadingColOrder[x])
  if(shouldConjugate)
  {
    tmp = 1L : nrow(mV)
    rst = lapply(rst, function(x) tmp[-x])
  }
  if(!fixedSize) rst = backout(rst, len)


  list(solution = rst, INT = c(INT, list(compressedDim = compressedDim)))
}













































