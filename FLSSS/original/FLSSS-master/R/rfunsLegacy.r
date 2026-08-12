

z_insert <- function(v, i, val)
{
  v = c(v, 0)
  # length(v) - 1L cannot be less than i because v's size now has been incremented
  v[i : length(v)] = c(val, v[i : (length(v) - 1L)])
  v
}


z_spread <- function(x) # x is number of target sums
{
  xrsv = 1L : x
  existx = c(1L, x, integer(x - 2L))
  for(i in 3L : length(existx))
  {
    remain = xrsv[-existx[1L : (i - 1L)]]
    existx[i] = remain[which.max(sapply(remain, function(u) min(abs(u - existx))))]
  }
  existx
}




z_reorderTarget <- function(keyTarget, givenOrder)
{
  if(length(keyTarget) == 1L) return(keyTarget)
  # reorder the keyTarget or not:
  if(length(givenOrder) != 1L)
  {
    if(length(givenOrder) < length(keyTarget))
    {
      return(paste0("randomizedTargetOrder should be at least of size ", length(keyTarget)))
    }
    keyTarget = keyTarget[givenOrder[givenOrder <= length(keyTarget)]]
  }
  else if(givenOrder)
  {
    keyTarget = sample(keyTarget, length(keyTarget))
  }
  else
  {
    keyTargetLen = length(keyTarget)
    # if(keyTargetLen %% 2L == 0L)
    # {
    #   givenOrder = as.integer(t(cbind((as.integer(keyTargetLen / 2L) + 1L) : keyTargetLen, as.integer(keyTargetLen / 2L) : 1L)))
    # }
    # else
    # {
    #   firstEle = as.integer(keyTargetLen / 2L) + 1L
    #   givenOrder = c(firstEle, as.integer(t(cbind((firstEle + 1L) : keyTargetLen, (firstEle - 1L) : 1L))))
    # }
    # keyTarget = keyTarget[givenOrder]
    keyTarget = keyTarget[z_spread(keyTargetLen)]
  }
  keyTarget
}




# The first column is the key column. The key column is always bounded from two sides
z_colOrder <- function(ncolV, dl, du) # dl and du describes things before the key column is added
{
  theOrder = 1L : ncolV
  dl = dl + 1L
  du = du + 1L
  keyPlace = ncolV - du + 1L
  # i = 1L
  if(keyPlace != 1L)
  {
    theOrder[1L : (keyPlace - 1L)] = theOrder[2L : keyPlace]
    theOrder[keyPlace] = 1L
  }
  list(theOrder = theOrder, zeroBasedKeyInd = keyPlace - 1L)
}




z_adhereIndexCol <- function(len, mV, ME, target, randomOrderKeyTarget, dl, du, verbose = F)
{
  # make key column
  tmpV = t(mV)
  key = integer(nrow(mV))
  i = 2L
  while(i <= length(key))
  {
    if(!all(tmpV[, i] >= tmpV[, i - 1L]))
    {
      key[i] = key[i - 1L] + 1L
    }
    else
    {
      key[i] = key[i - 1L]
    }
    i = i + 1L
  }
  rm(tmpV); gc()


  mV = data.frame(key, mV)
  k = unlist(lapply(mV[-1], function(x) max(0, max(-diff(x) / 1))))


  mV[-1] = mapply(function(x, y) x * mV[[1]] + y, k, mV[-1], SIMPLIFY = F)
  keyME = min(ME / k)
  if(is.nan(keyME)) keyME = 0
  # lsum = sum(0L : (len - 1L))
  # usum = sum((nrow(mV) - len) : (nrow(mV) - 1L))
  lsum = sum(mV[[1]][1L : len])
  usum = sum(mV[[1]][(nrow(mV) - len + 1L) : nrow(mV)])


  if(2 * keyME < 1) By = 1L
  else if(2L * as.integer(keyME) ==  as.integer(2 * keyME)) By = as.integer(2 * keyME) + 1L
  else By = as.integer(2 * keyME)


  keyTarget = seq(lsum + as.integer(keyME), usum, by = By)
  if(keyTarget[length(keyTarget)] < usum) keyTarget = c(keyTarget, usum)


  ME = ME - k * as.integer(keyME)
  if(as.integer(keyME) < 1L) keyME = 1e-1
  ME = c(keyME, ME)


  k = c(1, k)
  target = c(0, target)


  # Sort columns of mV, the scale factor, the target such that the key column is bounded from both sides
  tmpOrder = z_colOrder(ncol(mV), dl, du)
  mV = mV[tmpOrder$theOrder]
  k = k[tmpOrder$theOrder]
  target = target[tmpOrder$theOrder]
  ME = ME[tmpOrder$theOrder]


  if(verbose) cat("With index key column, atom tasks = ", length(keyTarget), "\n")


  reorderedTarget = z_reorderTarget(keyTarget, randomOrderKeyTarget)
  if(!is.integer(reorderedTarget) & !is.numeric(reorderedTarget)) return(reorderedTarget)
  list(v = mV, keyTarget = reorderedTarget, scaleFactor = k, originalTarget = target, ME = ME, dldu = c(dl + 1L, du + 1L), zeroBasedKeyInd = tmpOrder$zeroBasedKeyInd)
}




z_viaLeadingCol <- function(len, mV, ME, target, randomOrderKeyTarget, dl, du, verbose = F)
{
  intervalsNeeded = integer(ncol(mV))
  alphaList = list()
  for(i in 1L : ncol(mV))
  {
    v = mV[order(mV[[i]]), ]
    alpha = unlist(lapply(v[-i], function(x) max(0, max(-diff(x) / diff(v[[i]])))))
    alpha[is.nan(alpha)] = 0 # alpha is of size d - 1
    MEiub = ME[-i] / alpha
    tmpNintervals = as.integer(ME[i] / MEiub) + 1L
    while(!all(ME[i] / tmpNintervals <= MEiub)) tmpNintervals = tmpNintervals + 1L
    intervalsNeeded[i] = max(tmpNintervals)
    alphaList[[i]] = alpha
  }
  i = which.min(intervalsNeeded)
  minIntervals = as.integer(intervalsNeeded[i])


  if(verbose) cat("Without index key column, atom tasks = ", minIntervals, "\n")


  if(minIntervals > len * (nrow(mV) - len) + 1L) return(NULL)


  me = 2 * ME[i] / minIntervals / 2
  allKeyTarget = seq(target[i] - ME[i] + me, target[i] + (ME[i] - me), len = minIntervals)
  scaleFactor = alphaList[[i]]
  mV[-i] = mapply(function(x, y)
  {
    x + mV[[i]] * y
  }, mV[-i], scaleFactor, SIMPLIFY = F)
  ME[-i] = ME[-i] - scaleFactor * me
  ME[i] = me
  leadingColOrder = order(mV[[i]])
  mV = mV[leadingColOrder, ]


  # if(i != 1L)
  # {
  #   mV[1L : i] = mV[c(i, 1L : (i - 1L))]
  #   ME[1L : i] = ME[c(i, 1L : (i - 1L))]
  # } # Now mV's first column is the key column


  scaleFactor = z_insert(scaleFactor, i, 1) # dummy scaleFactor at position i
  target[i] = 0;
  dl = max(i, dl)
  du = max(ncol(mV) - i + 1L, du) # suboptimal method but effective, do not change it since there might be totally
  # useless columns in mV


  reorderedTarget = z_reorderTarget(allKeyTarget, randomOrderKeyTarget)
  if(!is.integer(reorderedTarget) & !is.numeric(reorderedTarget)) return(reorderedTarget)
  # list(v = mV, keyTarget = reorderedTarget, nonKeyOriginalTarget = target[-i], scaleFactor = scaleFactor, ME = ME, leadingColOrder = leadingColOrder)
  list(v = mV, keyTarget = reorderedTarget, originalTarget = target, scaleFactor = scaleFactor, ME = ME, leadingColOrder = leadingColOrder, dldu = c(dl, du), zeroBasedKeyInd = i - 1L)
}





































