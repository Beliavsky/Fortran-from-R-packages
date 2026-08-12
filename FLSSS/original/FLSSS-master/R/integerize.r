

{
## we find that scale and shift + integerize is equivalent to integerize + scale and shift
# shiftAndScale <- function(V)
# {
#   N = nrow(V)
#   d = ncol(V)
#   sft = pmax(0 - apply(v, 2, function(x) min(diff(x))), rep(0, d))
#   for(i in 1L : d)
#   {
#     V[, i] = V[, i] + sft[i] * (0L : (N - 1L))
#   }
#   V
# }


# z_comoIntegerV <- function(len, integerizeRst, outputColMajor = T)
# {
#   v = integerizeRst$integerized
#   target = integerizeRst$target
#
#
#   deltas = apply(v, 2, function(x)
#   {
#     tmp = min(diff(x))
#     if(tmp >= 0L) 0L
#     else 0L - tmp
#   })
#   tv = t(v)
#   mf = integer(ncol(tv))
#   for(i in 2L : ncol(tv))
#   {
#     if(!all(tv[, i] >= tv[, i - 1L])) mf[i] = mf[i - 1L] + 1L
#     else mf[i] = mf[i - 1L]
#   }
#   rm(tv); gc()
#   comoV = v
#   comoVtarget = list()
#   mfseq = seq(sum(mf[1L : len]), sum(mf[(nrow(v) - len + 1L) : nrow(v)]), by = 1L)
#   for(i in 1L : ncol(v))
#   {
#     comoV[, i] = comoV[, i] + deltas[i] * mf
#     comoVtarget[[i]] = mfseq * deltas[i] + target[i]
#   }
#   comoVtarget = as.matrix(as.data.frame(comoVtarget))
#   dimnames(comoVtarget) = NULL
#   tmp = comoV[1, ]
#   comoV = apply(comoV, 2, function(x) x - x[1])
#   comoVtarget = t(apply(comoVtarget, 1, function(x) x - tmp * len))
#   largestSubsetSum = apply(comoV, 2, function(x) sum(x[(length(x) - len + 1L) : length(x)]))
#   largestSubsetSum = pmax(largestSubsetSum, comoVtarget[nrow(comoVtarget), ])
#   Ndigits = as.integer(log2(largestSubsetSum)) + 1L + 1L
#   ME = integerizeRst$ME
#   {
#     low = sum(comoV[1L : len, 1])
#     high = sum(comoV[(nrow(comoV) - len + 1L) : nrow(comoV), 1])
#     r = (comoVtarget[1] - low) / (high - low)
#     tmpOrder = order(abs(1L : nrow(comoVtarget) - ((nrow(comoVtarget) - 1) * r + 1)))
#     comoVtarget = comoVtarget[tmpOrder, ]
#   }
#   comoV = comoV
#   if(outputColMajor) list(comoV = comoV, comoVtarget = comoVtarget, ME = integerizeRst$ME, Ndigits = Ndigits)
#   else list(comoV = t(comoV), comoVtarget = t(comoVtarget), ME = integerizeRst$ME, Ndigits = Ndigits)
# }
#
#
# z_integerizeComonotonize <- function(len, v, target, ME, outputColMajor = T)
# {
#   len = as.integer(len)
#   minV = apply(v, 2, function(x) min(x))
#   v = apply(v, 2, function(x) x - min(x))
#   target = target - minV * len
#   rst = z_integerize(len, v, target, ME)
#   minV = apply(rst$integerized, 2, function(x) min(x))
#   rst$integerized = apply(rst$integerized, 2, function(x) x - min(x))
#   rst$target = rst$target - minV * len
#   z_comoIntegerV(len, rst, outputColMajor)
# }
}




# # V is a matrix. ncol(V) is the dimensions
# z_integerize64bit <- function(len, V, target, ME)
# {
#   N = nrow(V)
#   largestSubsetSum = apply(V, 2, function(x) sum(sort(x)[(N - len + 1L) : N]))
#   target = pmin(largestSubsetSum, target)
#   tmp = z_integerize(len, V, target, ME)
#   V = tmp$integerized
#   target = tmp$target
#   ME = tmp$ME
#   integerized = list(V = V, target = target, ME = ME);
#   largestSubsetSum = apply(V, 2, function(x) sum(sort(x)[(N - len + 1L) : N]))
#   tmp = z_which64intAndSize(largestSubsetSum)
#   which64int = tmp$which64int
#   bitSize = tmp$bitSize
#   V = z_collapseTo64int(t(V), which64int, bitSize)
#   target = z_collapseTo64int(as.matrix(target), which64int, bitSize)
#   ME = z_collapseTo64int(as.matrix(ME), which64int, bitSize)
#   mask = z_mask(which64int, bitSize)
#   list(V = V, target = target, ME = ME
#        , which64int = which64int, bitSize = bitSize
#        , mask = mask
#        , integerized = integerized
#        )
# }




# mV, targetMat and mME are integers
z_crunchIntegers <- function(len, mV, targetMat, mME, dlst = NULL, dl = NULL, dust = NULL, du = NULL, maxMag = NULL)
{
  N = nrow(mV)
  if(is.null(maxMag)) maxMag = apply(mV, 2L, function(x) sum(sort(x)[(N - len + 1L) : N]))
  tmp = z_which64intAndSize(maxMag)
  which64int = tmp$which64int


  if(!is.null(dlst))
  {
    dlstnew = which64int[dlst + 1L]
    dlnew = which64int[dlst + dl - 1L + 1L] - dlstnew + 1L
    dustnew = which64int[dust + 1L]
    dunew = which64int[dust + du - 1L + 1L] - dustnew + 1L
  }


  bitSize = tmp$bitSize
  # print(bitSize)
  rst = list()
  rst$mV = t(z_collapseTo64int(t(mV), which64int, bitSize))
  rst$targetMat = z_collapseTo64int(as.matrix(targetMat), which64int, bitSize)
  rst$mME = as.numeric(z_collapseTo64int(as.matrix(mME), which64int, bitSize))
  rst$maskV = z_mask(which64int, bitSize)
  if(!is.null(dlst)) rst$dim = c(dlstnew, dlnew, dustnew, dunew)
  rst
}







































