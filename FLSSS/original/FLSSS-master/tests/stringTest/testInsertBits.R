

Rcpp::sourceCpp("src/arbitraryDimFLSSStests/insertBits.cpp", verbose = T)


while(T){


Nbuffer = 5L
maxN64bitIntPerbuffer = 5L
N32bitIntOfAllbuffers = lapply(1:Nbuffer, function(x)
  sample(maxN64bitIntPerbuffer, 1) * 2L)
X = lapply(N32bitIntOfAllbuffers, function(x)
{
  sample(2147483647L, x)
})
NeffectiveBits = unlist(lapply(X, function(x)
{
  sample(length(x) * 32L, 1)
}))
correctRst = as.integer(unlist(mapply(function(x, y)
{
  intToBits(x)[1:y]
}, X, NeffectiveBits, SIMPLIFY = F)))


rst = testInsertBits(X, NeffectiveBits)
rst = as.integer(intToBits(rst))
err = sum(abs(rst[1:length(correctRst)] - correctRst))
if(length(correctRst) < length(rst))
  err = err + sum(rst[(length(correctRst) + 1L):length(rst)])
if(err != 0) break
}
































