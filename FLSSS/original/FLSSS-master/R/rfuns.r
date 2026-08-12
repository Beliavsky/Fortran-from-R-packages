

z_seq <- function(x, y)
{
  if(x > y) integer(0)
  else x : y
}




# V is a column-major matrix and sorted by the leading column
z_mTargetMat <- function(len, mV, mTarget, mME, keyColInd, dl, du)
{
  N = nrow(mV)
  d = ncol(mV)


  tmp = t(mV)
  keyIndex = integer(N)
  for(i in 2L : N)
  {
    if(!all(tmp[, i - 1] <= tmp[, i])) keyIndex[i] = keyIndex[i - 1] + 1L
    else keyIndex[i] = keyIndex[i - 1]
  }
  rm(tmp); gc()
  # print(keyIndex)


  # attach index column
  {
    # the index column would be placed at dl + 1L
    if(keyColInd > dl)
    {
      keyColInd = keyColInd + 1L
    }
    mV = cbind(mV[, 1L : dl], numeric(N), mV[, z_seq(dl + 1L, d)])
    mME = c(mME[1L : dl], 1 / 8, mME[z_seq(dl + 1L, d)])
    mTarget = c(mTarget[1L : dl], 0, mTarget[z_seq(dl + 1L, d)])
    dl = dl + 1L
    du = du + 1L
    d = d + 1L
  }


  scaleFactor = apply(mV, 2, function(x)
  {
    tmp = min(diff(x))
    if(tmp >= 0) 0
    else 0 - tmp
  })
  scaleFactor[dl] = 1


  mV = mV + as.matrix(as.data.frame(lapply(scaleFactor, function(x) keyIndex * x)))
  keyTargets = sum(keyIndex[1L : len]) : sum(keyIndex[(N - len + 1L) : N])
  targetMat = as.matrix(as.data.frame(lapply(keyTargets, function(x) scaleFactor * x + mTarget)))
  {
    lowS = sum(mV[1L : len, keyColInd])
    highS = sum(mV[(N - len + 1L) : N, keyColInd])
    targetPercentile = (mTarget[keyColInd] - lowS) / (highS - lowS)
    pivot = (keyTargets[length(keyTargets)] - keyTargets[1]) * targetPercentile + keyTargets[1]
  }
  targetMat = targetMat[, order(abs(keyTargets - pivot))]
  dimnames(targetMat) = NULL
  list(mV = mV, targetMat = targetMat, mME = mME, d = d, dl = dl, du = du, keyColInd = keyColInd)
}




z_mTargetMatNoKey <- function(len, mV, mTarget, mME, ordering, dl, du)
{
  N = nrow(mV)
  d = ncol(mV)


  tmp = t(mV)
  keyIndex = integer(N)
  for(i in 2L : N)
  {
    if(!all(tmp[, i - 1] <= tmp[, i])) keyIndex[i] = keyIndex[i - 1] + 1L
    else keyIndex[i] = keyIndex[i - 1]
  }
  rm(tmp); gc()


  # attach index column
  {
    mV = cbind(mV[, 1L : dl], numeric(N), mV[, z_seq(dl + 1L, d)])
    mME = c(mME[1L : dl], 1 / 8, mME[z_seq(dl + 1L, d)])
    mTarget = c(mTarget[1L : dl], 0, mTarget[z_seq(dl + 1L, d)])
    dl = dl + 1L
    du = du + 1L
    d = d + 1L
  }


  scaleFactor = apply(mV, 2, function(x)
  {
    tmp = min(diff(x))
    if(tmp >= 0) 0
    else 0 - tmp
  })
  scaleFactor[dl] = 1


  mV = mV + as.matrix(as.data.frame(lapply(scaleFactor, function(x) keyIndex * x)))
  keyTargets = sum(keyIndex[1L : len]) : sum(keyIndex[(N - len + 1L) : N])
  targetMat = as.matrix(as.data.frame(lapply(keyTargets, function(x) scaleFactor * x + mTarget)))
  if(!is.null(ordering))
  {
    ordering = unique(ordering %% N + 1L)
    targetMat = targetMat[, ordering]
  }
  else
  {
    pivot = (keyTargets[1] + keyTargets[length(keyTargets)]) / 2
    targetMat = targetMat[, order(abs(keyTargets - pivot))]
  }
  dimnames(targetMat) = NULL
  list(mV = mV, targetMat = targetMat, mME = mME, dl = dl, du = du, d = d)
}




z_mTargetMatForKnapsack <- function(len, mV, mTarget, mME)
{
  N = nrow(mV)
  tmp = t(mV)
  keyIndex = integer(N)
  for(i in 2L : N)
  {
    if(!all(tmp[, i - 1] <= tmp[, i])) keyIndex[i] = keyIndex[i - 1] + 1L
    else keyIndex[i] = keyIndex[i - 1]
  }
  rm(tmp); gc()


  # attach index column
  {
    mV = cbind(numeric(N), mV)
    mME = c(1 / 8, mME)
    mTarget = c(0, mTarget)
  }


  scaleFactor = apply(mV, 2, function(x)
  {
    tmp = min(diff(x))
    if(tmp >= 0) 0
    else 0 - tmp
  })
  scaleFactor[1] = 1


  mV = mV + as.matrix(as.data.frame(lapply(scaleFactor, function(x) keyIndex * x)))
  keyTargets = sum(keyIndex[(N - len + 1L) : N]) : sum(keyIndex[1L : len])
  targetMat = as.matrix(as.data.frame(lapply(keyTargets, function(x) scaleFactor * x + mTarget)))
  dimnames(targetMat) = NULL
  list(mV = mV, targetMat = targetMat, mME = mME)
}








#--------------------------------------------------------------------------------------------------


# V is a column-major matrix and sorted by the leading column
z_mTargetMatINT <- function(len, mV, mTarget, mME, keyColInd, dl, du)
{
  N = nrow(mV)
  d = ncol(mV)


  tmp = t(mV)
  keyIndex = integer(N)
  for(i in 2L : N)
  {
    if(!all(tmp[, i - 1L] <= tmp[, i])) keyIndex[i] = keyIndex[i - 1L] + 1L
    else keyIndex[i] = keyIndex[i - 1L]
  }
  rm(tmp); gc()


  # attach index column
  {
    # the index column would be placed at dl + 1L
    if(keyColInd > dl)
    {
      keyColInd = keyColInd + 1L
    }
    mV = cbind(mV[, 1L : dl], integer(N), mV[, z_seq(dl + 1L, d)])
    mME = c(mME[1L : dl], 0L, mME[z_seq(dl + 1L, d)])
    mTarget = c(mTarget[1L : dl], 0L, mTarget[z_seq(dl + 1L, d)])
    dl = dl + 1L
    du = du + 1L
    d = d + 1L
  }


  scaleFactor = apply(mV, 2L, function(x)
  {
    tmp = min(diff(x))
    if(tmp >= 0L) 0L
    else 0L - tmp
  })
  scaleFactor[dl] = 1L


  mV = mV + as.matrix(as.data.frame(lapply(scaleFactor, function(x) keyIndex * x)))
  keyTargets = sum(keyIndex[1L : len]) : sum(keyIndex[(N - len + 1L) : N])
  targetMat = as.matrix(as.data.frame(lapply(keyTargets, function(x) scaleFactor * x + mTarget)))


  # reorder targetMat. Put the target that is the most likely to yield a solution first
  {
    lowS = sum(mV[1L : len, keyColInd])
    highS = sum(mV[(N - len + 1L) : N, keyColInd])
    targetPercentile = (mTarget[keyColInd] - lowS) / (highS - lowS)
    pivot = (keyTargets[length(keyTargets)] - keyTargets[1L]) * targetPercentile + keyTargets[1L]
  }
  targetMat = targetMat[, order(abs(keyTargets - pivot))]
  dimnames(targetMat) = NULL
  list(mV = mV, targetMat = targetMat, mME = mME, d = d, dl = dl, du = du, keyColInd = keyColInd)
}




z_mTargetMatNoKeyINT <- function(len, mV, mTarget, mME, ordering, dl, du)
{
  N = nrow(mV)
  d = ncol(mV)


  tmp = t(mV)
  keyIndex = integer(N)
  for(i in 2L : N)
  {
    if(!all(tmp[, i - 1L] <= tmp[, i])) keyIndex[i] = keyIndex[i - 1L] + 1L
    else keyIndex[i] = keyIndex[i - 1L]
  }
  rm(tmp); gc()


  # attach index column
  {
    mV = cbind(mV[, 1L : dl], integer(N), mV[, z_seq(dl + 1L, d)])
    mME = c(mME[1L : dl], 0L, mME[z_seq(dl + 1L, d)])
    mTarget = c(mTarget[1L : dl], 0L, mTarget[z_seq(dl + 1L, d)])
    dl = dl + 1L
    du = du + 1L
    d = d + 1L
  }


  scaleFactor = apply(mV, 2L, function(x)
  {
    tmp = min(diff(x))
    if(tmp >= 0L) 0L
    else 0L - tmp
  })
  scaleFactor[dl] = 1L


  mV = mV + as.matrix(as.data.frame(lapply(scaleFactor, function(x) keyIndex * x)))
  keyTargets = sum(keyIndex[1L : len]) : sum(keyIndex[(N - len + 1L) : N])
  targetMat = as.matrix(as.data.frame(lapply(keyTargets, function(x) scaleFactor * x + mTarget)))
  if(!is.null(ordering))
  {
    ordering = unique(ordering %% N + 1L)
    targetMat = targetMat[, ordering]
  }
  else
  {
    pivot = (keyTargets[1] + keyTargets[length(keyTargets)]) / 2L
    targetMat = targetMat[, order(abs(keyTargets - pivot))]
  }
  dimnames(targetMat) = NULL
  list(mV = mV, targetMat = targetMat, mME = mME, dl = dl, du = du, d = d)
}




z_mTargetMatForKnapsackINT <- function(len, mV, mTarget, mME)
{
  N = nrow(mV)


  tmp = t(mV)
  keyIndex = integer(N)
  for(i in 2L : N)
  {
    if(!all(tmp[, i - 1L] <= tmp[, i])) keyIndex[i] = keyIndex[i - 1L] + 1L
    else keyIndex[i] = keyIndex[i - 1L]
  }
  rm(tmp); gc()


  # attach index column
  {
    mV = cbind(integer(N), mV)
    mME = c(0L, mME)
    mTarget = c(0L, mTarget)
  }


  scaleFactor = apply(mV, 2L, function(x)
  {
    tmp = min(diff(x))
    if(tmp >= 0L) 0L
    else 0L - tmp
  })
  scaleFactor[1L] = 1L


  mV = mV + as.matrix(as.data.frame(lapply(scaleFactor, function(x) keyIndex * x)))
  keyTargets = sum(keyIndex[(N - len + 1L) : N]) : sum(keyIndex[1L : len])
  targetMat = as.matrix(as.data.frame(lapply(keyTargets, function(x) scaleFactor * x + mTarget)))
  dimnames(targetMat) = NULL
  list(mV = mV, targetMat = targetMat, mME = mME)
}








# mV has been comonotonized and ordered
# each column of targetMat is a target
z_filterTargetFindLargestMagnitude <- function(len, mV, targetMat, mME, lowestSS = NULL, largestSS = NULL)
{
  N = nrow(mV)
  u = (N - len + 1L) : N
  if(is.null(largestSS)) largestSS = apply(mV, 2, function(x) sum(x[u]))
  u = 1L : len
  if(is.null(lowestSS)) lowestSS = apply(mV, 2, function(x) sum(x[u]))
  # targetMat[, i] - mME should not be greater than largestSS
  # targetMat[, i] + mME should not be less than lowestSS
  L = lowestSS - mME
  U = largestSS + mME
  targetMat = as.matrix(targetMat[, apply(targetMat, 2, function(x) all(x <= U & x >= L))]) # 'as.matrix()' here ensures 'targetMat' stays as a matrix, not a fucking vector
  maxMag = pmax(apply(targetMat, 1, function(x) max(abs(x))) + mME, abs(L), abs(U))
  list(maxMag = maxMag, targetMat = targetMat)
}




z_filterByBoundsOneDim <- function(len, V, target, LB, UB)
{
  od = 1L : length(V)
  VtouchedInd = unique(unlist(apply(cbind(LB, UB), 1, function(x) x[1] : x[2])))
  Vtouched = V[VtouchedInd]
  originalIndTouched = od[VtouchedInd]
  LUB = apply(cbind(LB, UB), 1, function(x)
  {
    match(c(x[1], x[2]), originalIndTouched)
  })
  LB = LUB[1, ]
  UB = LUB[2, ]
  mapFun = function(solution, mapIndex) {mapIndex[solution]}
  list(V = Vtouched, LB = LB, UB = UB, mapIndex = originalIndTouched, mapFun = mapFun)
}







































